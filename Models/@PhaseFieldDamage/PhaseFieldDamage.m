classdef PhaseFieldDamage < BaseModel
    %Resolves the phase field evolution, and mechanical momentum balances.
	%inputs required:
	% physics_in{1}.type = "PhaseFieldDamage";
    % physics_in{1}.Egroup = "Internal";
    % physics_in{1}.young = 200e9;	%Youngs Modulus [Pa]
    % physics_in{1}.poisson = 0.3;	%Poisson ratio [-]
	% physics_in{1}.kmin = 1e-10;		%residual stiffness factor [-]
	% physics_in{1}.l=l;				%Phase-field length scale [m]
	% physics_in{1}.Gc=2.0e3;			%Fracture release energy [J/m^2]
	% physics_in{1}.GDegrade = 0.9;	%Maximum hydrogen degradation factor
	% physics_in{1}.NL = 1e6;			%Concentration of interstitial lattice sites [mol/m^3]
	% physics_in{1}.gb = 30e3;		%Grain boundary binding energy [J/mol]
    
    properties
        mesh			%Pointer to mesh object
        myName			%String with the name of this model
        myGroup			%String indicating the element group this model operates on
        myGroupIndex	%Index of element group
        dofSpace		%Pointer to degree of freedom object
        dofTypeIndices	%Indices of degrees of freedom associated with this model
        
        poisson	%Poisson ratio [-]
        young	%youngs modulus [Pa]
        D_el	%Linear-elastic-plane-strain-stiffness matrix
		kmin	%residual stiffness [-]
		l		%phase-field length scale [m]
		Gc		%Fracture release energy
		GDegrade%fracture energy degradation rate due to hydrogen
		NL		%concentration of interstitial lattice sites
		gb		%binding energy of grain boundaries

		dx_Step		%step in which displacements are resolved
		phi_step	%step in which phasefield variable is resolved
		CL_step		%step in which lattice hydrogen concentration is resolved

		doInit		%Flag to indicate whether first-time initialization has been performed

		Hist		%Eleastic energy history parameter of current increment
		HistOld		%Eleastic energy history parameter of previous increment

		LFrac		%Estimated crack length

		T = 293.15;
		R = 8.3144;
    end
    
    methods
        function obj = PhaseFieldDamage(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "PhaseFieldDamage";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            [obj.dofTypeIndices, stp]  = obj.dofSpace.getDofType({"dx","dy","phi","CL"});
			obj.dx_Step = stp(1);
			obj.phi_step = stp(3);
			obj.CL_step = stp(4);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            %% constuct plane-strain elastic stiffness matrix
            obj.poisson = inputs.poisson;
            obj.young = inputs.young;
			obj.kmin = inputs.kmin;
			obj.l = inputs.l;
			obj.Gc = inputs.Gc;
			obj.GDegrade = inputs.GDegrade;
			obj.NL = inputs.NL;
			obj.gb = inputs.gb;
            
            D_el = zeros(4,4);
            a = obj.young / ((1.0 + obj.poisson) * (1.0 - 2.0*obj.poisson));

            D_el(1, 1) = a * (1.0 - obj.poisson);
            D_el(1, 2) = a * obj.poisson;
            D_el(1, 3) = a * obj.poisson;
            D_el(2, 1) = a * obj.poisson;
            D_el(2, 2) = a * (1.0 - obj.poisson);
            D_el(2, 3) = a * obj.poisson;
            D_el(3, 1) = a * obj.poisson;
            D_el(3, 2) = a * obj.poisson;
            D_el(3, 3) = a * (1.0 - obj.poisson);
            D_el(4, 4) = a * 0.5 * (1.0 - 2.0*obj.poisson);
            
            obj.D_el = D_el;

			obj.Hist = zeros(size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1), obj.mesh.ipcount1D^2);
			obj.HistOld = obj.Hist;

			obj.doInit = true;
		end

		function OncePerStep(obj, physics, stp)
			if (stp == obj.phi_step && obj.doInit)
				Hdom = max(obj.mesh.Nodes(:,2));
				Lfrac = 5e-3;
	
				%% set Values
				nodecons = [];

            	allNodes = obj.mesh.GetAllNodesForGroup(obj.myGroupIndex);
				PhiDofs = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), allNodes);
            	initvals = 0*allNodes;
				for i=1:length(allNodes)
					xy = [obj.mesh.Nodes(allNodes(i),1), obj.mesh.Nodes(allNodes(i),2)];
					
					if (xy(1)<Lfrac)
						dst = Hdom/2-xy(2);
					else
						dst=sqrt((xy(2)-Hdom/2)^2+(xy(1)-Lfrac)^2);
					end
					
					pf = exp(-abs(dst)/(obj.l));  

					initvals(i) = pf;
				end

				physics.StateVec{obj.phi_step}(PhiDofs) = initvals;
				physics.StateVec_Old{obj.phi_step}(PhiDofs) = initvals;

				%% set history field
				for n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                	Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                	[N, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
					G2 = obj.mesh.getG2(obj.myGroupIndex, n_el);

					dofsPhi = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
					PHI = physics.StateVec{obj.phi_step}(dofsPhi);

					for ip=1:length(w)
						NPhi = N(ip,:)*PHI;
						GPhi = squeeze(G(ip,:,:))'*PHI;

						d_dam_fun = -2*(1-obj.kmin)*(1-NPhi);
						H = abs((NPhi/(obj.l)+obj.l*(GPhi'*GPhi))/(d_dam_fun+obj.kmin));
						
						obj.Hist(n_el, ip) = H;
						obj.HistOld(n_el, ip) = H;
					end		

				end

				obj.doInit = false;
			end

		end
        
        function Commit(obj, physics, commit_type)
			if (commit_type == "Pathdep")
				obj.HistOld = obj.Hist;
			end
        end

        function getKf(obj, physics, stp)

			if (stp == obj.dx_Step)
            	fprintf("        PhaseFieldDamage get Matrix:")
            	t = tic;
            	
            	dofmatX = [];
            	dofmatY = [];
            	kmat = [];

				fvec = [];
            	dofvec = [];

				SVec = physics.StateVec;

            	parfor n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                	Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                	[N, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);

                	dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                	dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                	dofsXY = [dofsX; dofsY];
					dofsPhi = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
					dofsCL = obj.dofSpace.getDofIndices(obj.dofTypeIndices(4), Elem_Nodes);

                	X = SVec{obj.dx_Step}(dofsX);
                	Y = SVec{obj.dx_Step}(dofsY);
                	XY = [X;Y];
					PHI = SVec{obj.phi_step}(dofsPhi);
					CL = SVec{obj.CL_step}(dofsCL);

                	f_el = zeros(length(dofsXY), 1);
                	K_el = zeros(length(dofsXY));
                	for ip=1:length(w)
                    	B = obj.getB(G(ip,:,:));
                    	strain = B*XY;

						ff = min(max(N(ip,:)*PHI,0),1);
						dam_fun = (1-obj.kmin)*(1-ff)^2+obj.kmin;

                    	stress = dam_fun*obj.D_el*strain;
                    	f_el = f_el + B'*stress*w(ip);
                    	K_el = K_el + B'*dam_fun*obj.D_el*B*w(ip);
                	end

                	[dofmatxloc,dofmatyloc] = ndgrid(dofsXY,dofsXY);
                	dofmatX = [dofmatX; dofmatxloc(:)];
                	dofmatY = [dofmatY; dofmatyloc(:)];
                	kmat = [kmat; K_el(:)];

					fvec = [fvec; f_el];
                	dofvec = [dofvec; dofsXY];
            	end 

				physics.fint{stp} = physics.fint{stp} + sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint{stp}), 1);
            	physics.K{stp} = physics.K{stp} + sparse(dofmatX, dofmatY, kmat, length(physics.fint{stp}),length(physics.fint{stp}));
            	
            	tElapsed = toc(t);
            	fprintf("            (Assemble time:"+string(tElapsed)+")\n");
			end

			if (stp == obj.phi_step)
            	fprintf("        PhaseFieldDamage get Matrix:")
            	t = tic;
            	
            	dofmatX = [];
            	dofmatY = [];
            	kmat = [];

				SVec = physics.StateVec;

				fvec = [];
            	dofvec = [];

				maxIP = size(obj.Hist,2);
				HNew = zeros(size(obj.Hist,1),size(obj.Hist,2));
				
				fraccheck = 0;
            	parfor n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)     %parfor
                	Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                	[N, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);

                	dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                	dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                	dofsXY = [dofsX; dofsY];
					dofsPhi = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
					dofsCL = obj.dofSpace.getDofIndices(obj.dofTypeIndices(4), Elem_Nodes);

                	X = SVec{obj.dx_Step}(dofsX);
                	Y = SVec{obj.dx_Step}(dofsY);
                	XY = [X;Y];
					PHI = SVec{obj.phi_step}(dofsPhi);
					CL = SVec{obj.CL_step}(dofsCL);

                	f_el = zeros(length(dofsPhi), 1);
                	K_el = zeros(length(dofsPhi));
                	for ip=1:maxIP
						gamma = 1/(2*obj.l)*(N(ip,:)*PHI)^2+obj.l/2*((squeeze(G(ip,:,:))'*PHI)'*(squeeze(G(ip,:,:))'*PHI));

						% resistance
						f_el = f_el + w(ip)*( N(ip,:)'*N(ip,:)*PHI/(obj.l) ...
											 + obj.l*squeeze(G(ip,:,:))*squeeze(G(ip,:,:))'*PHI);

						K_el = K_el + w(ip)*( N(ip,:)'*N(ip,:)/(obj.l) ...
											 + obj.l*squeeze(G(ip,:,:))*squeeze(G(ip,:,:))');

						% driving force
                    	B = obj.getB(G(ip,:,:));
                    	strain = B*XY;
						Gc_loc = obj.Gc * ( 1-obj.GDegrade*max(0,(N(ip,:)*CL/obj.NL)/(N(ip,:)*CL/obj.NL+exp(-obj.gb/obj.T/obj.R))) );
						Energy_el = 0.5*strain'*obj.D_el*strain/Gc_loc;

						dam_fun = (1-obj.kmin)*(1-N(ip,:)*PHI)^2+obj.kmin;
						d_dam_fun = -2*(1-obj.kmin)*(1-N(ip,:)*PHI);
						d_dam_dphi = 2*(1-obj.kmin);
						
						H = obj.HistOld(n_el, ip);
						if (H>=Energy_el)
							f_el = f_el + d_dam_fun*w(ip)*N(ip,:)'*H;

							K_el = K_el + d_dam_dphi*w(ip)*H*N(ip,:)'*N(ip,:);

							HNew(n_el, ip) = obj.HistOld(n_el, ip);
						else
							f_el = f_el + d_dam_fun*w(ip)*N(ip,:)'*Energy_el;

							K_el = K_el + w(ip)*d_dam_dphi*Energy_el*N(ip,:)'*N(ip,:);

							HNew(n_el, ip) = Energy_el;
						end

						fraccheck = fraccheck + gamma*w(ip);
                	end

                	[dofmatxloc,dofmatyloc] = ndgrid(dofsPhi,dofsPhi);
                	dofmatX = [dofmatX; dofmatxloc(:)];
                	dofmatY = [dofmatY; dofmatyloc(:)];
                	kmat = [kmat; K_el(:)];

					fvec = [fvec; f_el];
                	dofvec = [dofvec; dofsPhi];
            	end 

				physics.fint{stp} = physics.fint{stp} + sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint{stp}), 1);
            	physics.K{stp} = physics.K{stp} + sparse(dofmatX, dofmatY, kmat, length(physics.fint{stp}),length(physics.fint{stp}));
            	
				obj.LFrac = fraccheck;

            	tElapsed = toc(t);
            	fprintf("            (Assemble time:"+string(tElapsed)+")\n");
				fprintf("Fracture Length: "+string(fraccheck)+"\n");

				obj.Hist = HNew;
			end
        end
        
        function B = getB(~, grads)
            cp_count = size(grads, 2);
            B = zeros(4, cp_count*2);
            for ii = 1:cp_count %using plane strain e_zz = 0
				%dx
				B(1, ii) = grads(1,ii, 1);
				B(4, ii) = grads(1,ii, 2);

				%dy
				B(2, ii + cp_count) = grads(1,ii, 2);
				B(4, ii + cp_count) = grads(1,ii, 1);
            end
        end

        function [hasInfo, provided] = Provide_Info(obj, physics, var, elems, loc)
           hasInfo = false;
           provided = [];
           
           if (var == "stresses" || var=="sxx" || var=="syy" || var=="szz" || var=="sxy" || var=="sh") %elem, ip, componenets
               hasInfo = true;
               if (var == "stresses")
                   provided = zeros(length(elems), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount, 4);
               else
                   provided = zeros(length(elems), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount);
               end
                for el=1:length(elems)
                    Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, elems(el));
                    [N, G, w] = obj.mesh.getVals(obj.myGroupIndex, elems(el));

                    dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                    dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
					dofsPhi = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
					
                    X = physics.StateVec{obj.dx_Step}(dofsX);
                    Y = physics.StateVec{obj.dx_Step}(dofsY);
                    XY = [X;Y];
					PHI = physics.StateVec{obj.phi_step}(dofsPhi);

                    for ip=1:length(w)
						ff = min(max(N(ip,:)*PHI,0),1);
						dam_fun = (1-obj.kmin)*(1-ff)^2+obj.kmin;

                        B = obj.getB(G(ip,:,:));
                        strain = B*XY;
                        stress = dam_fun*obj.D_el*strain;

                        if (loc == "Interior")
                            switch var 
                                case "stresses"
                                    provided(el, ip, :) = stress;
                                case "sxx"
                                    provided(el, ip) = stress(1);
                                case "syy"
                                    provided(el, ip) = stress(2);
                                case "szz"
                                    provided(el, ip) = stress(3);
                                case "sxy"
                                    provided(el, ip) = stress(4);
								case "sh"
                                    provided(el, ip) = (stress(1)+stress(2)+stress(3))/3;
                            end
                        end
                        
                    end
                end
           end
        end
        
    end
end

