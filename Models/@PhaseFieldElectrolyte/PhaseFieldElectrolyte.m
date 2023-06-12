classdef PhaseFieldElectrolyte < BaseModel
    %ElectrolyteFrac Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
		fracTypeIndices

        D
		z
		k
		Lumped
		pH0
		NaCl

		F_const = 96485.3329;
		R_const = 8.31446261815324;
		T_const = 293.15;

		n_species

		C_step
		dx_step
		phi_step

		NAds
		ksurf
		NL
		Em
		Lumpedsurf
		l

		phi_trigger
		h0
		Flowtype  %  "WuLorenzis" / "subgrid"
		fracCons

		Heights
		Normals

		ActiveElems
		fracNodes
    end
    
    methods
        function obj = PhaseFieldElectrolyte(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "PhaseFieldElectrolyte";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            [obj.dofTypeIndices, stp] = obj.dofSpace.getDofType({"dx","dy","phi","CL"});
			obj.dx_step = stp(1);
			obj.phi_step = stp(3);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));

			[obj.fracTypeIndices, stp] = obj.dofSpace.getDofType({"Epot", "theta", "H", "OH","Na", "Cl","Fe","FeOH"});
			obj.C_step = stp(7);

            
            %% get parameters
            obj.D = inputs.D;
			obj.pH0 = inputs.pH0;
			obj.NaCl = inputs.NaCl;
			obj.z = inputs.z;
			obj.k=inputs.k;
			obj.Lumped=inputs.Lumped;

			obj.NAds = inputs.NAds;
			obj.ksurf = inputs.ksurf;
			obj.NL = inputs.NL;
			obj.Em = inputs.Em;
			obj.Lumpedsurf = inputs.Lumpedsurf;
			obj.h0 = inputs.h0;
			obj.Flowtype = inputs.Flowtype;
			obj.phi_trigger = inputs.phi_trigger;

			obj.n_species = length(obj.D);
			obj.l = inputs.l;

			obj.ActiveElems = [];
			obj.fracNodes = [];
        end
        
		function OncePerStep(obj, physics, stp)
			if (stp == obj.C_step)
				%Detect active elements
				NewNodes = [];
				NewElems = [];

				for n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
					Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
					dofsPhi = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
					Phi = physics.StateVec{obj.phi_step}(dofsPhi);
					if (max(abs(Phi))>obj.phi_trigger)
						if isempty(find(n_el==obj.ActiveElems, 1))
							NewElems = [NewElems, n_el];
							obj.ActiveElems = [obj.ActiveElems, n_el];
					
							for i=1:length(Elem_Nodes)
								if isempty(find(Elem_Nodes(i)==obj.fracNodes, 1))
									NewNodes = [NewNodes, Elem_Nodes(i)];
									obj.fracNodes = [obj.fracNodes, Elem_Nodes(i)];
								end
	
							end
						end
	
					end
				end
	
				%add relevant fracture specific dofs
				obj.dofSpace.addDofs(obj.fracTypeIndices, NewNodes);
	
				dofsE = obj.dofSpace.getDofIndices(obj.fracTypeIndices(1), NewNodes);
				dofsT = obj.dofSpace.getDofIndices(obj.fracTypeIndices(2), NewNodes);
				for i=1:obj.n_species
					dofsC(:,i) = obj.dofSpace.getDofIndices(obj.fracTypeIndices(i+2), NewNodes);
				end
	
				initH = 1000*10^(-obj.pH0);
				initOH = 1000*10^(-14+obj.pH0);
				initCl = obj.NaCl;
				initNa = initCl-initH+initOH;
	
				initState = [initH, initOH, initNa, initCl, 0, 0];
	
				physics.StateVec{obj.C_step}(dofsE) = 0;	physics.StateVec_Old{obj.C_step}(dofsE) = 0;
				physics.StateVec{obj.C_step}(dofsT) = 0;	physics.StateVec_Old{obj.C_step}(dofsT) = 0;
				for i=1:obj.n_species
					physics.StateVec{obj.C_step}(dofsC(:,i)) = initState(i);	
					physics.StateVec_Old{obj.C_step}(dofsC(:,i)) = initState(i);
				end
	
				%construct constrain vectors
				clear dofsE dofsT dofsC

				dofsE = obj.dofSpace.getDofIndices(obj.fracTypeIndices(1), obj.fracNodes);
				dofsT = obj.dofSpace.getDofIndices(obj.fracTypeIndices(2), obj.fracNodes);
				for i=1:obj.n_species
					dofsC(:,i) = obj.dofSpace.getDofIndices(obj.fracTypeIndices(i+2), obj.fracNodes);
				end

				obj.fracCons = zeros(0,2);
				for i=1:length(obj.fracNodes)
					xy = [obj.mesh.Nodes(obj.fracNodes(i),1), obj.mesh.Nodes(obj.fracNodes(i),2)];
					if (xy(1)==0)
						addcon = [dofsE(i), 0];
						obj.fracCons = [obj.fracCons; addcon];
	
						for j=1:obj.n_species
							addcon = [dofsC(i,j), initState(j)];
							obj.fracCons = [obj.fracCons; addcon];
						end
					end
				end

				%% opening height map
				obj.Construct_H_map(physics);
			end
		end

		function Construct_H_map(obj, physics)
			fprintf("        PhaseFieldElectrolyte Mapping Heights")

			obj.Heights = zeros(length(obj.ActiveElems), obj.mesh.ipcount1D^2);
			obj.Normals = zeros(length(obj.ActiveElems), obj.mesh.ipcount1D^2,2);

			xy = zeros(length(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems), obj.mesh.ipcount1D^2, 2);
			u = zeros(length(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems), obj.mesh.ipcount1D^2, 2);
			dphi = zeros(length(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems), obj.mesh.ipcount1D^2, 2);

			Svec = physics.StateVec;
			xmax = 0; ymax = 0;
			for el=1:length(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems)
				[N, G, ~] = obj.mesh.getVals(obj.myGroupIndex, el);
				coords = obj.mesh.getIPCoords(obj.myGroupIndex, el);

        		Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, el);
        		dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
        		dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
				dofsPhi = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);

				X = Svec{obj.dx_step}(dofsX);
                Y = Svec{obj.dx_step}(dofsY);
				PHI= Svec{obj.phi_step}(dofsPhi);

				for ip=1:obj.mesh.ipcount1D^2
					xy(el,ip,1:2) = coords(:,ip);
					u(el,ip,1:2)  = [N(ip,:)*X N(ip,:)*Y];
					dphi(el,ip,:) = squeeze(G(ip,:,:))'*PHI;
				end

				xmax = max(max(coords(1,:)),xmax);
				ymax = max(max(coords(2,:)),ymax);
			end

			%% assuming solely horizontal propagation
			x = xy(:,:,1); y = xy(:,:,2);
			ux = u(:,:,1); uy = u(:,:,2);
			dphix = squeeze(dphi(:,:,1)); dphiy = squeeze(dphi(:,:,2));

			F_ux = scatteredInterpolant(x(:),y(:),ux(:));
			 warning('off','last')
			F_uy = scatteredInterpolant(x(:),y(:),uy(:));
			F_fx = scatteredInterpolant(x(:),y(:),dphix(:));
			F_fy = scatteredInterpolant(x(:),y(:),dphiy(:));

			x_eval = linspace(0,xmax,1000);  ny = 1000;
			h_est = zeros(size(x_eval));
			for iy=0:ny
				y_eval = 0*x_eval+iy*ymax./ny;
				h_est = h_est + abs(F_ux(x_eval,y_eval).*F_fx(x_eval,y_eval)*0 + F_uy(x_eval,y_eval).*F_fy(x_eval,y_eval))*ymax./ny;
			end
			h_est = h_est;
			Fh = griddedInterpolant(x_eval, h_est);
			
			for n_el=1:length(obj.ActiveElems)
				el = obj.ActiveElems(n_el);
				coords = obj.mesh.getIPCoords(obj.myGroupIndex, el);
				
				for ip=1:obj.mesh.ipcount1D^2
					coords_ip = coords(:,ip);
					obj.Heights(n_el,ip) = Fh(coords_ip(1));
					obj.Normals(n_el,ip,1) = F_fx(coords_ip(1),coords_ip(2));
					obj.Normals(n_el,ip,2) = F_fy(coords_ip(1),coords_ip(2));
					obj.Normals(n_el,ip,:) = obj.Normals(n_el,ip,:)/(sqrt(obj.Normals(n_el,ip,1)^2+obj.Normals(n_el,ip,2)^2));

					if obj.Heights(n_el,ip)<1e-8
						obj.Heights(n_el,ip) = 1e-8;
					end
				end
			end

			obj.plotHeightData(x_eval, h_est);

			%% still to actually calculate
% 			obj.Heights = zeros(length(obj.ActiveElems), obj.mesh.ipcount1D^2)+obj.h0;
% 			obj.Normals = zeros(length(obj.ActiveElems), obj.mesh.ipcount1D^2,2);
% 			obj.Normals(:,:,1) = 0;
% 			obj.Normals(:,:,2) = 1;
		end

        function getKf(obj, physics, stp)
			if (stp == obj.C_step)
            	fprintf("        PhaseFieldElectrolyte get Matrix:")
            	t = tic;
            	
            	dt = physics.dt;
	
				%constrain required dofs
            	physics.condofs{obj.C_step} = [physics.condofs{obj.C_step}; obj.fracCons(:,1)];
            	physics.convals{obj.C_step} = [physics.convals{obj.C_step}; obj.fracCons(:,2)];

				%stiffness matrix
            	dofmatX = [];
            	dofmatY = [];
            	kmat = [];
            	fvec = [];
            	dofvec = [];
	
				fraccheck = 0;

				Svec = physics.StateVec;
            	SvecOld = physics.StateVec_Old;


            	parfor n_el=1:length(obj.ActiveElems)   %PARFOR
					el = obj.ActiveElems(n_el);

                	Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, el);
                	[N, G, w] = obj.mesh.getVals(obj.myGroupIndex, el);
	
					dofsE = obj.dofSpace.getDofIndices(obj.fracTypeIndices(1), Elem_Nodes);
					dofsC = zeros(length(dofsE),obj.n_species);
					for i=1:obj.n_species
						dofsC(:,i) = obj.dofSpace.getDofIndices(obj.fracTypeIndices(i+2), Elem_Nodes);
					end
					dofsT = obj.dofSpace.getDofIndices(obj.fracTypeIndices(2), Elem_Nodes);

					dofsPHI= obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
					dofsCL= obj.dofSpace.getDofIndices(obj.dofTypeIndices(4), Elem_Nodes);
	
					C = Svec{obj.C_step}(dofsC);
					E = Svec{obj.C_step}(dofsE);
					T = Svec{obj.C_step}(dofsT);

					PHI = Svec{obj.phi_step}(dofsPHI);
					CL = Svec{obj.C_step}(dofsCL);
	
					C_OLD = SvecOld{obj.C_step}(dofsC);
					E_OLD = SvecOld{obj.C_step}(dofsE);
					T_OLD = SvecOld{obj.C_step}(dofsT);

					PHI_OLD = SvecOld{obj.phi_step}(dofsPHI);
					CL_OLD = SvecOld{obj.C_step}(dofsCL);

					q_C    = zeros(length(dofsE), obj.n_species);
					dqC_dC = zeros(length(dofsE), length(dofsE), obj.n_species, obj.n_species);
					dqC_dE = zeros(length(dofsE), length(dofsE), obj.n_species);
					dqC_dT = zeros(length(dofsE), length(dofsE), obj.n_species);

					q_E = zeros(length(dofsE), 1);	
					dqE_dE = zeros(length(dofsE), length(dofsE));
					dqE_dC = zeros(length(dofsE), length(dofsE), obj.n_species);

					q_T = zeros(length(dofsE), 1);
					dqT_dT = zeros(length(dofsE), length(dofsE));
					dqT_dE = zeros(length(dofsE), length(dofsE));
					dqT_dC = zeros(length(dofsE), length(dofsE), obj.n_species);
					dqT_dCL = zeros(length(dofsE), length(dofsCL));

					q_CL = zeros(length(dofsCL), 1);
					dqCL_dCL = zeros(length(dofsCL), length(dofsCL));
					dqCL_dT = zeros(length(dofsCL), length(dofsE));

	
					C_Lumped_Vol = zeros(length(dofsE), 1);
					C_Lumped_surf = zeros(length(dofsE), 1);
					for ip=1:length(w)

						%phase-field distributor function
						Nphi = N(ip,:)*PHI;
						Gphi = squeeze(G(ip,:,:))'*PHI;
						gamma = max(0,1/(2*obj.l)*Nphi^2+obj.l/2*(Gphi'*Gphi));
	
						distributor_D = [0,0];
						distributor_D_offset = [0,0];
						if (obj.Flowtype == "WuLorenzis") %macro-scale description
							distributor_cap = Nphi;
							distributor_cap_offset = Nphi+obj.h0;
							m=1;
							distributor_D(1)= Nphi^m;
							distributor_D(2)= 0;%Nphi^m;
							distributor_D_offset(1)= Nphi^m+obj.h0;
							distributor_D_offset(2)= obj.h0;%Nphi^m+obj.h0;
						else
							distributor_cap = gamma*obj.Heights(n_el, ip);
							distributor_cap_offset = gamma*obj.Heights(n_el, ip)+obj.h0;
							distributor_D(1)= gamma*obj.Heights(n_el, ip);
							distributor_D(2)= gamma*1000;
							distributor_D_offset(1) = gamma*obj.Heights(n_el, ip)+obj.h0;
							distributor_D_offset(2)= gamma*1000+obj.h0;
						end

						%capacity
						for s=1:obj.n_species
							q_C(:,s) = q_C(:,s) + distributor_cap*w(ip)/dt*N(ip,:)'*N(ip,:)*(C(:,s)-C_OLD(:,s));
							dqC_dC(:,:,s,s) = dqC_dC(:,:,s,s) + distributor_cap_offset*w(ip)/dt*N(ip,:)'*N(ip,:);
						end

						%diffusion
						Df = zeros(2,2);
						DfNum = zeros(2,2);
						for s=1:obj.n_species
							Df(1,1) = distributor_D(1)*obj.D(s);
							Df(2,2) = distributor_D(2)*obj.D(s);
							DfNum(1,1) = distributor_D_offset(1)*obj.D(s);
							DfNum(2,2) = distributor_D_offset(2)*obj.D(s);
	
							q_C(:,s) = q_C(:,s) + w(ip)*squeeze(G(ip,:,:))*Df*squeeze(G(ip,:,:))'*C(:,s);
							dqC_dC(:,:,s,s) = dqC_dC(:,:,s,s) + w(ip)*squeeze(G(ip,:,:))*DfNum*squeeze(G(ip,:,:))';
						end

						%electroneutrality
						Df = zeros(2,2);
						for s=1:obj.n_species
							Df(1,1) = distributor_D(1)*obj.D(s);
							Df(2,2) = distributor_D(2)*obj.D(s);

							DfNum(1,1) = distributor_D_offset(1)*obj.D(s);
							DfNum(2,2) = distributor_D_offset(2)*obj.D(s);
	
							fct = obj.z(s)*obj.F_const/obj.R_const/obj.T_const;
							
							q_C(:,s) =	q_C(:,s) + w(ip)*(N(ip,:)*C(:,s))*squeeze(G(ip,:,:))*(fct*Df*squeeze(G(ip,:,:))'*E);		
							dqC_dC(:,:,s,s) = dqC_dC(:,:,s,s) + w(ip)*squeeze(G(ip,:,:))*(fct*DfNum*squeeze(G(ip,:,:))'*E)*N(ip,:);
							dqC_dE(:,:,s) = dqC_dE(:,:,s) + w(ip)*(N(ip,:)*C(:,s))*squeeze(G(ip,:,:))*fct*DfNum*squeeze(G(ip,:,:))';
	
							q_E = q_E + distributor_cap*w(ip)*obj.z(s)*N(ip,:)'*N(ip,:)*C(:,s);
							dqE_dC(:,:,s) = dqE_dC(:,:,s) + distributor_cap_offset*w(ip)*obj.z(s)*N(ip,:)'*N(ip,:);
						end

						% surface capacity
						q_T = q_T + w(ip)*2*gamma*obj.NAds/dt*N(ip,:)'*N(ip,:)*(T-T_OLD);
						dqT_dT = dqT_dT+ w(ip)*2*gamma*obj.NAds/dt*N(ip,:)'*N(ip,:);

						%% volume reactions
	
						% water
						if (obj.Lumped(1)==false)
							rwc = obj.k(1);
							rw = rwc*(1e-8-(N(ip,:)*C(:,1))*(N(ip,:)*C(:,2)));
							drw_dH  = rwc*(-(N(ip,:))*(N(ip,:)*C(:,2)));
							drw_dOH = rwc*(-(N(ip,:)*C(:,1))*(N(ip,:)));

							q_C(:,1) = q_C(:,1) - distributor_cap*w(ip)*N(ip,:)'*rw;
							q_C(:,2) = q_C(:,2) - distributor_cap*w(ip)*N(ip,:)'*rw;
		
							dqC_dC(:,:,1,1) = dqC_dC(:,:,1,1) - distributor_cap*w(ip)*N(ip,:)'*drw_dH;
							dqC_dC(:,:,1,2) = dqC_dC(:,:,1,2) - distributor_cap*w(ip)*N(ip,:)'*drw_dOH;
							dqC_dC(:,:,2,1) = dqC_dC(:,:,2,1) - distributor_cap*w(ip)*N(ip,:)'*drw_dH;
							dqC_dC(:,:,2,2) = dqC_dC(:,:,2,2) - distributor_cap*w(ip)*N(ip,:)'*drw_dOH;
						end
	
						%Fe
						if (obj.Lumped(2)==false)
							rcc = [obj.k(2), obj.k(3), obj.k(4)];
							rc = rcc(1)*(N(ip,:)*C(:,5))-rcc(2)*(N(ip,:)*C(:,6))*(N(ip,:)*C(:,1));
							drc_dH = -rcc(2)*(N(ip,:)*C(:,6))*N(ip,:);
							drc_dFE = rcc(1)*N(ip,:);
							drc_dFEOH = -rcc(2)*(N(ip,:)*C(:,1))*N(ip,:);
		
							rc2 = rcc(3)*(N(ip,:)*C(:,6));
							drc2_dFEOH = rcc(3)*N(ip,:);
	
							q_C(:,1) = q_C(:,1) - distributor_cap*w(ip)*N(ip,:)'*(rc+rc2);
							q_C(:,5) = q_C(:,5) - distributor_cap*w(ip)*N(ip,:)'*-rc;
							q_C(:,6) = q_C(:,6) - distributor_cap*w(ip)*N(ip,:)'*(rc-rc2);
		
							dqC_dC(:,:,1,1) = dqC_dC(:,:,1,1) - distributor_cap*w(ip)*N(ip,:)'*drc_dH;
							dqC_dC(:,:,1,5) = dqC_dC(:,:,1,5) - distributor_cap*w(ip)*N(ip,:)'*drc_dFE;
							dqC_dC(:,:,1,6) = dqC_dC(:,:,1,6) - distributor_cap*w(ip)*N(ip,:)'*(drc_dFEOH+drc2_dFEOH);
		
							dqC_dC(:,:,5,1) = dqC_dC(:,:,5,1) - distributor_cap*w(ip)*N(ip,:)'*-drc_dH;
							dqC_dC(:,:,5,5) = dqC_dC(:,:,5,5) - distributor_cap*w(ip)*N(ip,:)'*-drc_dFE;
							dqC_dC(:,:,5,6) = dqC_dC(:,:,5,6) - distributor_cap*w(ip)*N(ip,:)'*-drc_dFEOH;
		
							dqC_dC(:,:,6,1) = dqC_dC(:,:,6,1) - distributor_cap*w(ip)*N(ip,:)'*drc_dH;
							dqC_dC(:,:,6,5) = dqC_dC(:,:,6,5) - distributor_cap*w(ip)*N(ip,:)'*drc_dFE;
							dqC_dC(:,:,6,6) = dqC_dC(:,:,6,6) - distributor_cap*w(ip)*N(ip,:)'*(drc_dFEOH-drc2_dFEOH);
						end
	
						%% surface reactions
						CH = N(ip,:)*C(:,1);
						COH = N(ip,:)*C(:,2);
						CFE = N(ip,:)*C(:,5);
						theta = N(ip,:)*T;
						phil = N(ip,:)*E;
						CLat = N(ip,:)*CL;
	
						[react, dreact, products] = obj.reactions(CH, COH, CFE, theta, phil, CLat);
	
						% total reactions
						for r=1:7
							q_C(:,1) = q_C(:,1) - 2*gamma*w(ip)*N(ip,:)'*(react(r,1)-react(r,2))*products(r,1) *(1-obj.Lumpedsurf(r));
							q_C(:,2) = q_C(:,2) - 2*gamma*w(ip)*N(ip,:)'*(react(r,1)-react(r,2))*products(r,2)*(1-obj.Lumpedsurf(r));
							q_C(:,5) = q_C(:,5) - 2*gamma*w(ip)*N(ip,:)'*(react(r,1)-react(r,2))*products(r,3)*(1-obj.Lumpedsurf(r));
							q_T      = q_T      - 2*gamma*w(ip)*N(ip,:)'*(react(r,1)-react(r,2))*products(r,4) *(1-obj.Lumpedsurf(r));
							q_CL     = q_CL     - 2*gamma*w(ip)*N(ip,:)'*(react(r,1)-react(r,2))*products(r,5)*(1-obj.Lumpedsurf(r));
	
							for n=1:obj.n_species
								dqC_dC(:,:,1,n) = dqC_dC(:,:,1,n)-2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,1)*(1-obj.Lumpedsurf(r));
								dqC_dC(:,:,2,n) = dqC_dC(:,:,2,n)-2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,2)*(1-obj.Lumpedsurf(r));
								dqC_dC(:,:,5,n) = dqC_dC(:,:,5,n)-2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,3)*(1-obj.Lumpedsurf(r));
								dqT_dC(:,:,n)   = dqT_dC(:,:,n)  -2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,4)*(1-obj.Lumpedsurf(r));
							end
							dqC_dE(:,:,1)   = dqC_dE(:,:,1) - 2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,1)-dreact(r,2,1))*products(r,1)*(1-obj.Lumpedsurf(r));
							dqC_dE(:,:,2)   = dqC_dE(:,:,2) - 2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,1)-dreact(r,2,1))*products(r,2)*(1-obj.Lumpedsurf(r));
							dqC_dE(:,:,5)   = dqC_dE(:,:,5) - 2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,1)-dreact(r,2,1))*products(r,3)*(1-obj.Lumpedsurf(r));
							dqT_dE          = dqT_dE        - 2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,1)-dreact(r,2,1))*products(r,4)*(1-obj.Lumpedsurf(r));
	
							dqC_dT(:,:,1)   = dqC_dT(:,:,1) - 2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,2)-dreact(r,2,2))*products(r,1)*(1-obj.Lumpedsurf(r));
							dqC_dT(:,:,2)   = dqC_dT(:,:,2) - 2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,2)-dreact(r,2,2))*products(r,2)*(1-obj.Lumpedsurf(r));
							dqC_dT(:,:,5)   = dqC_dT(:,:,5) - 2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,2)-dreact(r,2,2))*products(r,3)*(1-obj.Lumpedsurf(r));
							dqT_dT          = dqT_dT        - 2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,2)-dreact(r,2,2))*products(r,4)*(1-obj.Lumpedsurf(r));
							dqCL_dT         = dqCL_dT       - 2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,2)-dreact(r,2,2))*products(r,5)*(1-obj.Lumpedsurf(r));
	
							dqT_dCL         = dqT_dCL       - 2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,3)-dreact(r,2,3))*products(r,4)*(1-obj.Lumpedsurf(r));
							dqCL_dCL        = dqCL_dCL      - 2*gamma*w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,3)-dreact(r,2,3))*products(r,5)*(1-obj.Lumpedsurf(r));
						end

						%% lumped integration factor
						C_Lumped_Vol = C_Lumped_Vol + distributor_cap*w(ip)*N(ip,:)';
						C_Lumped_surf= C_Lumped_surf+ gamma*w(ip)*N(ip,:)';
						fraccheck = fraccheck + gamma*w(ip);
					end
	
	
					for i=1:length(dofsE) %% lumped integrations			
						if (obj.Lumped(1)) %water
							rwc = obj.k(1);
							kdum = 1e10;
							phlim=1e-16;
	
							rw = rwc*(1e-8-(max(phlim,C(i,1)))*(max(phlim,C(i,2))));
							drw_dH  = rwc*(-(max(phlim,C(i,2))));
							drw_dOH = rwc*(-(max(phlim,C(i,1))));
	
							if (C(i,1)<phlim)
								rw = rw + (phlim*1.1-C(i,1))*kdum;
								drw_dH = drw_dH - kdum;
							end
							if (C(i,2)<phlim)
								rw = rw + (phlim*1.1-C(i,2))*kdum;
								drw_dOH = drw_dOH - kdum;
							end
	
							q_C(i,1) = q_C(i,1) - C_Lumped_Vol(i)*rw;
							q_C(i,2) = q_C(i,2) - C_Lumped_Vol(i)*rw;
		
							dqC_dC(i,i,1,1) = dqC_dC(i,i,1,1) - C_Lumped_Vol(i)*drw_dH;
							dqC_dC(i,i,1,2) = dqC_dC(i,i,1,2) - C_Lumped_Vol(i)*drw_dOH;
							dqC_dC(i,i,2,1) = dqC_dC(i,i,2,1) - C_Lumped_Vol(i)*drw_dH;
							dqC_dC(i,i,2,2) = dqC_dC(i,i,2,2) - C_Lumped_Vol(i)*drw_dOH;
						end
						if (obj.Lumped(2)) %Fe
							rcc = [obj.k(2), obj.k(3), obj.k(4)];
							rc = rcc(1)*C(i,5)-rcc(2)*C(i,6)*C(i,1);
							drc_dH = -rcc(2)*C(i,6);
							drc_dFE = rcc(1);
							drc_dFEOH = -rcc(2)*C(i,1);
		
							rc2 = rcc(3)*C(i,6);
							drc2_dFEOH = rcc(3);
	
							q_C(i,1) = q_C(i,1) - C_Lumped_Vol(i)*(rc+rc2);
							q_C(i,5) = q_C(i,5) - C_Lumped_Vol(i)*-rc;
							q_C(i,6) = q_C(i,6) - C_Lumped_Vol(i)*(rc-rc2);
		
							dqC_dC(i,i,1,1) = dqC_dC(i,i,1,1) - C_Lumped_Vol(i)*drc_dH;
							dqC_dC(i,i,1,5) = dqC_dC(i,i,1,5) - C_Lumped_Vol(i)*drc_dFE;
							dqC_dC(i,i,1,6) = dqC_dC(i,i,1,6) - C_Lumped_Vol(i)*(drc_dFEOH+drc2_dFEOH);
		
							dqC_dC(i,i,5,1) = dqC_dC(i,i,5,1) - C_Lumped_Vol(i)*-drc_dH;
							dqC_dC(i,i,5,5) = dqC_dC(i,i,5,5) - C_Lumped_Vol(i)*-drc_dFE;
							dqC_dC(i,i,5,6) = dqC_dC(i,i,5,6) - C_Lumped_Vol(i)*-drc_dFEOH;
		
							dqC_dC(i,i,6,1) = dqC_dC(i,i,6,1) - C_Lumped_Vol(i)*drc_dH;
							dqC_dC(i,i,6,5) = dqC_dC(i,i,6,5) - C_Lumped_Vol(i)*drc_dFE;
							dqC_dC(i,i,6,6) = dqC_dC(i,i,6,6) - C_Lumped_Vol(i)*(drc_dFEOH-drc2_dFEOH);
						end


						[react, dreact, products] = obj.reactions(C(i,1), C(i,2),C(i,5), T(i), E(i), CL(i));
	
						for r=1:7
							q_C(i,1) = q_C(i,1) - 2*C_Lumped_surf(i)*(react(r,1)-react(r,2))*products(r,1)*obj.Lumpedsurf(r);
							q_C(i,2) = q_C(i,2) - 2*C_Lumped_surf(i)*(react(r,1)-react(r,2))*products(r,2)*obj.Lumpedsurf(r);
							q_C(i,5) = q_C(i,5) - 2*C_Lumped_surf(i)*(react(r,1)-react(r,2))*products(r,3)*obj.Lumpedsurf(r);
							q_T(i)   = q_T(i)   - 2*C_Lumped_surf(i)*(react(r,1)-react(r,2))*products(r,4)*obj.Lumpedsurf(r);
							q_CL(i)  = q_CL(i)  - 2*C_Lumped_surf(i)*(react(r,1)-react(r,2))*products(r,5)*obj.Lumpedsurf(r);
	
							for n=1:obj.n_species
								dqC_dC(i,i,1,n) = dqC_dC(i,i,1,n)-2*C_Lumped_surf(i)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,1)*obj.Lumpedsurf(r);
								dqC_dC(i,i,2,n) = dqC_dC(i,i,2,n)-2*C_Lumped_surf(i)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,2)*obj.Lumpedsurf(r);
								dqC_dC(i,i,5,n) = dqC_dC(i,i,5,n)-2*C_Lumped_surf(i)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,3)*obj.Lumpedsurf(r);
								dqT_dC(i,i,n)   = dqT_dC(i,i,n)  -2*C_Lumped_surf(i)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,4)*obj.Lumpedsurf(r);
							end
							dqC_dE(i,i,1)   = dqC_dE(i,i,1) - 2*C_Lumped_surf(i)*(dreact(r,1,1)-dreact(r,2,1))*products(r,1)*obj.Lumpedsurf(r);
							dqC_dE(i,i,2)   = dqC_dE(i,i,2) - 2*C_Lumped_surf(i)*(dreact(r,1,1)-dreact(r,2,1))*products(r,2)*obj.Lumpedsurf(r);
							dqC_dE(i,i,5)   = dqC_dE(i,i,5) - 2*C_Lumped_surf(i)*(dreact(r,1,1)-dreact(r,2,1))*products(r,3)*obj.Lumpedsurf(r);
							dqT_dE(i,i)     = dqT_dE(i,i)   - 2*C_Lumped_surf(i)*(dreact(r,1,1)-dreact(r,2,1))*products(r,4)*obj.Lumpedsurf(r);
	
							dqC_dT(i,i,1)   = dqC_dT(i,i,1) - 2*C_Lumped_surf(i)*(dreact(r,1,2)-dreact(r,2,2))*products(r,1)*obj.Lumpedsurf(r);
							dqC_dT(i,i,2)   = dqC_dT(i,i,2) - 2*C_Lumped_surf(i)*(dreact(r,1,2)-dreact(r,2,2))*products(r,2)*obj.Lumpedsurf(r);
							dqC_dT(i,i,5)   = dqC_dT(i,i,5) - 2*C_Lumped_surf(i)*(dreact(r,1,2)-dreact(r,2,2))*products(r,3)*obj.Lumpedsurf(r);
							dqT_dT(i,i)     = dqT_dT(i,i)   - 2*C_Lumped_surf(i)*(dreact(r,1,2)-dreact(r,2,2))*products(r,4)*obj.Lumpedsurf(r);
							dqCL_dT(i,i)    = dqCL_dT(i,i)  - 2*C_Lumped_surf(i)*(dreact(r,1,2)-dreact(r,2,2))*products(r,5)*obj.Lumpedsurf(r);
	
							dqT_dCL(i,i)    = dqT_dCL(i,i)  - 2*C_Lumped_surf(i)*(dreact(r,1,3)-dreact(r,2,3))*products(r,4)*obj.Lumpedsurf(r);
							dqCL_dCL(i,i)   = dqCL_dCL(i,i) - 2*C_Lumped_surf(i)*(dreact(r,1,3)-dreact(r,2,3))*products(r,5)*obj.Lumpedsurf(r);
						end

					end
	
					for s1=1:obj.n_species
						for s2=1:obj.n_species
							[dofmatxloc,dofmatyloc] = ndgrid(dofsC(:,s1),dofsC(:,s2));
							dofmatX = [dofmatX; dofmatxloc(:)];
							dofmatY = [dofmatY; dofmatyloc(:)];
							tmp = dqC_dC(:,:,s1,s2);
							kmat = [kmat; tmp(:)];
						end
	
						[dofmatxloc,dofmatyloc] = ndgrid(dofsC(:,s1),dofsE);
						dofmatX = [dofmatX; dofmatxloc(:)];
						dofmatY = [dofmatY; dofmatyloc(:)];
						tmp = squeeze(dqC_dE(:,:,s1));
						kmat = [kmat; tmp(:)];
	
						[dofmatxloc,dofmatyloc] = ndgrid(dofsE,dofsC(:,s1));
						dofmatX = [dofmatX; dofmatxloc(:)];
						dofmatY = [dofmatY; dofmatyloc(:)];
						tmp = squeeze(dqE_dC(:,:,s1));
						kmat = [kmat; tmp(:)];

						[dofmatxloc,dofmatyloc] = ndgrid(dofsC(:,s1),dofsT);
						dofmatX = [dofmatX; dofmatxloc(:)];
						dofmatY = [dofmatY; dofmatyloc(:)];
						tmp = dqC_dT(:,:,s1);
						kmat = [kmat; tmp(:)];

						[dofmatxloc,dofmatyloc] = ndgrid(dofsT,dofsC(:,s1));
						dofmatX = [dofmatX; dofmatxloc(:)];
						dofmatY = [dofmatY; dofmatyloc(:)];
						tmp = dqT_dC(:,:,s1);
						kmat = [kmat; tmp(:)];
	
						fvec = [fvec; q_C(:,s1)];
						dofvec = [dofvec; dofsC(:,s1)];
					end
	
					[dofmatxloc,dofmatyloc] = ndgrid(dofsE,dofsE);
                	dofmatX = [dofmatX; dofmatxloc(:)];
                	dofmatY = [dofmatY; dofmatyloc(:)];
                	kmat = [kmat; dqE_dE(:)];
	
					[dofmatxloc,dofmatyloc] = ndgrid(dofsT,dofsT);
                	dofmatX = [dofmatX; dofmatxloc(:)];
                	dofmatY = [dofmatY; dofmatyloc(:)];
                	kmat = [kmat; dqT_dT(:)];
	
					[dofmatxloc,dofmatyloc] = ndgrid(dofsT,dofsE);
                	dofmatX = [dofmatX; dofmatxloc(:)];
                	dofmatY = [dofmatY; dofmatyloc(:)];
                	kmat = [kmat; dqT_dE(:)];
	
					[dofmatxloc,dofmatyloc] = ndgrid(dofsT,dofsCL);
                	dofmatX = [dofmatX; dofmatxloc(:)];
                	dofmatY = [dofmatY; dofmatyloc(:)];
                	kmat = [kmat; dqT_dCL(:)];
	
					[dofmatxloc,dofmatyloc] = ndgrid(dofsCL,dofsCL);
                	dofmatX = [dofmatX; dofmatxloc(:)];
                	dofmatY = [dofmatY; dofmatyloc(:)];
                	kmat = [kmat; dqCL_dCL(:)];
	
					[dofmatxloc,dofmatyloc] = ndgrid(dofsCL,dofsT);
                	dofmatX = [dofmatX; dofmatxloc(:)];
                	dofmatY = [dofmatY; dofmatyloc(:)];
                	kmat = [kmat; dqCL_dT(:)];

					fvec = [fvec; q_E];
                	dofvec = [dofvec; dofsE];

					fvec = [fvec; q_T];
					dofvec = [dofvec; dofsT];

					fvec = [fvec; q_CL];
					dofvec = [dofvec; dofsCL];
            	end 
				
            	physics.fint{obj.C_step} = physics.fint{obj.C_step} + sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint{obj.C_step}), 1);
            	physics.K{obj.C_step} = physics.K{obj.C_step} + sparse(dofmatX, dofmatY, kmat, length(physics.fint{obj.C_step}),length(physics.fint{obj.C_step}));
	
            	tElapsed = toc(t);
            	fprintf("            (Assemble time:"+string(tElapsed)+")\n");
				fprintf("		Fracture Length: "+string(fraccheck)+"\n");
			end
		end

        
		function [react, dreact, products] = reactions(obj, CH, COH, CFE, theta, phil, CLat)
					Cmax = 1e3;

					products(:,1) = [-1, -1, 0, 0, 0, 0, 0];    %H+
					products(:,2) = [0, 0, 0, 0, 1, 1, 0];      %OH-
					products(:,3) = [0, 0, 0, 0, 0, 0, -1];      %Fe
					products(:,4) = [1, -1, -2, -1, 1, -1, 0];  %theta
					products(:,5) = [0, 0, 0, 1, 0, 0, 0];      %C_L

					react = zeros(7,2); %reaction, forward/backward
					dreact= zeros(7,2,3+obj.n_species); %reaction, forward/backwards, E/T/C_L/species

					dCH = 1;
					if (CH<0)
						CH = 0;
						dCH = 0;
					end
					dCOH = 1;
					if (COH<0)
						COH = 0;
						dCOH = 0;
					end

% 					dCFE = 1;
% 					if (CFE>Cmax)
% 						CFE = Cmax;
% 						dCFE = 0;
% 					end

					% reaction 1 
					react(1,1)    = obj.ksurf(1,1)*CH *(1-theta)*exp(-obj.ksurf(1,3)*(obj.Em-phil-obj.ksurf(1,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(1,1,2) = obj.ksurf(1,1)*CH *-1              *exp(-obj.ksurf(1,3)*(obj.Em-phil-obj.ksurf(1,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(1,1,4) = obj.ksurf(1,1)*dCH*(1-theta)*exp(-obj.ksurf(1,3)*(obj.Em-phil-obj.ksurf(1,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(1,1,1) = obj.ksurf(1,1)*CH *(1-theta)*exp(-obj.ksurf(1,3)*(obj.Em-phil-obj.ksurf(1,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.ksurf(1,3)*(-1)*obj.F_const/obj.R_const/obj.T_const);

					react(1,2)    = obj.ksurf(1,2)*theta *exp((1-obj.ksurf(1,3))*(obj.Em-phil-obj.ksurf(1,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(1,2,2) = obj.ksurf(1,2)              *exp((1-obj.ksurf(1,3))*(obj.Em-phil-obj.ksurf(1,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(1,2,1) = obj.ksurf(1,2)*theta *exp((1-obj.ksurf(1,3))*(obj.Em-phil-obj.ksurf(1,4))*obj.F_const/obj.R_const/obj.T_const)*((1-obj.ksurf(1,3))*(-1)*obj.F_const/obj.R_const/obj.T_const);

					% reaction 2
					react(2,1)    = obj.ksurf(2,1)*CH *theta*exp(-obj.ksurf(2,3)*(obj.Em-phil-obj.ksurf(2,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(2,1,4) = obj.ksurf(2,1)*dCH*theta*exp(-obj.ksurf(2,3)*(obj.Em-phil-obj.ksurf(2,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(2,1,2) = obj.ksurf(2,1)*CH              *exp(-obj.ksurf(2,3)*(obj.Em-phil-obj.ksurf(2,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(2,1,1) = obj.ksurf(2,1)*CH *theta*exp(-obj.ksurf(2,3)*(obj.Em-phil-obj.ksurf(2,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.ksurf(2,3)*(-1)*obj.F_const/obj.R_const/obj.T_const);

					%react(2,2) = 0;

					%reaction 3
					react(3,1)    = obj.ksurf(3,1)*abs(theta)*theta;
					dreact(3,1,2) = 2*obj.ksurf(3,1)*abs(theta);

					%react(3,2) = 0;

					%reaction 4
					%if theta>=0
						react(4,1)    = obj.ksurf(4,1)*(obj.NL-max(0,CLat))*theta;
						dreact(4,1,2) = obj.ksurf(4,1)*(obj.NL-max(0,CLat));
						dreact(4,1,3) = obj.ksurf(4,1)*(-1)*theta;
					%end
					react(4,2)    =  obj.ksurf(4,2)*CLat*(1-theta);
					dreact(4,2,2) =  obj.ksurf(4,2)*CLat*(-1);
					dreact(4,2,3) =  obj.ksurf(4,2)*(1-theta);


					%reaction 5
					react(5,1)    = obj.ksurf(5,1)*(1-theta)*exp(-obj.ksurf(5,3)*(obj.Em-phil-obj.ksurf(5,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(5,1,2) = obj.ksurf(5,1)*-1       *exp(-obj.ksurf(5,3)*(obj.Em-phil-obj.ksurf(5,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(5,1,1) = obj.ksurf(5,1)*(1-theta)*exp(-obj.ksurf(5,3)*(obj.Em-phil-obj.ksurf(5,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.ksurf(5,3)*(-1)*obj.F_const/obj.R_const/obj.T_const);

					react(5,2)    = obj.ksurf(5,2)*COH *theta*exp((1-obj.ksurf(5,3))*(obj.Em-phil-obj.ksurf(5,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(5,2,2) = obj.ksurf(5,2)*COH       *exp((1-obj.ksurf(5,3))*(obj.Em-phil-obj.ksurf(5,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(5,2,1) = obj.ksurf(5,2)*COH *theta*exp((1-obj.ksurf(5,3))*(obj.Em-phil-obj.ksurf(5,4))*obj.F_const/obj.R_const/obj.T_const)*((1-obj.ksurf(5,3))*(-1)*obj.F_const/obj.R_const/obj.T_const);
					dreact(5,2,5) = obj.ksurf(5,2)*dCOH*theta*exp((1-obj.ksurf(5,3))*(obj.Em-phil-obj.ksurf(5,4))*obj.F_const/obj.R_const/obj.T_const);

					%reaction 6
					react(6,1)    = obj.ksurf(6,1)*theta*exp(-obj.ksurf(6,3)*(obj.Em-phil-obj.ksurf(6,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(6,1,2) = obj.ksurf(6,1)      *exp(-obj.ksurf(6,3)*(obj.Em-phil-obj.ksurf(6,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(6,1,1) = obj.ksurf(6,1)*theta*exp(-obj.ksurf(6,3)*(obj.Em-phil-obj.ksurf(6,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.ksurf(6,3)*(-1)*obj.F_const/obj.R_const/obj.T_const);

					%react(6,2) = 0

					%corrosion
					react(7,1)    = obj.ksurf(7,1)*CFE*exp(-obj.ksurf(7,3)*(obj.Em-phil-obj.ksurf(7,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(7,1,1) = obj.ksurf(7,1)*CFE*exp(-obj.ksurf(7,3)*(obj.Em-phil-obj.ksurf(7,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.ksurf(7,3)*(-1)*obj.F_const/obj.R_const/obj.T_const);
					dreact(7,1,8) = obj.ksurf(7,1)*exp(-obj.ksurf(7,3)*(obj.Em-phil-obj.ksurf(7,4))*obj.F_const/obj.R_const/obj.T_const);
			
					react(7,2)    = obj.ksurf(7,2)*exp((1-obj.ksurf(7,3))*(obj.Em-phil-obj.ksurf(7,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(7,2,1) = obj.ksurf(7,2)*exp((1-obj.ksurf(7,3))*(obj.Em-phil-obj.ksurf(7,4))*obj.F_const/obj.R_const/obj.T_const)*((1-obj.ksurf(7,3))*(-1)*obj.F_const/obj.R_const/obj.T_const);
			
		end








		function plotHeightData(obj, x_eval, h_est)
% 			obj.Heights = zeros(length(obj.ActiveElems), obj.mesh.ipcount1D^2)+obj.h0;
% 			obj.Normals = zeros(length(obj.ActiveElems), obj.mesh.ipcount1D^2,2);

			figure(43191)
			clf
			subplot(3,1,1)
			plot(x_eval*1000,h_est*1000);
			hold on
			xlabel('x [mm]')
			ylabel('h [mm]')

			for n_el=1:length(obj.ActiveElems)
				el = obj.ActiveElems(n_el);
				coords = obj.mesh.getIPCoords(obj.myGroupIndex, el);

				x(n_el,:) = coords(1,:);
                y(n_el,:) = coords(2,:);
			end
			xlims = [min(min(x)), max(max(x))];
            ylims = [min(min(y)), max(max(y))];

			subplot(3,1,2)
            [xi,yi] = meshgrid(linspace(xlims(1), xlims(2), 1000), linspace(ylims(1), ylims(2), 1000));
            zi = griddata(x,y,obj.Heights,xi,yi,'linear');
            %s = pcolor(xi,yi,zi);
			s = surf(xi, yi, zi);
            s.EdgeColor = 'none';
            colorbar

			subplot(3,1,3)
            [xi,yi] = meshgrid(linspace(xlims(1), xlims(2), 100), linspace(ylims(1), ylims(2), 100));
            zix = griddata(x,y,squeeze(obj.Normals(:,:,1)),xi,yi,'linear');
			ziy = griddata(x,y,squeeze(obj.Normals(:,:,2)),xi,yi,'linear');
			quiver(xi,yi,zix,ziy)

			drawnow();

		end












		function plotFields(obj, physics)
			for el=1:length(obj.ActiveElems)
                elnodes =physics.mesh.Elementgroups{obj.myGroupIndex}.Elems(obj.ActiveElems(el),:);
				Edofs = physics.dofSpace.getDofIndices(obj.fracTypeIndices(1), elnodes);
				for s=1:obj.n_species
					Cdofs(:,s) = physics.dofSpace.getDofIndices(obj.fracTypeIndices(2+s), elnodes);
				end
				dofsT = obj.dofSpace.getDofIndices(obj.fracTypeIndices(2), elnodes);
				dofsPHI= obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), elnodes);
				dofsCL= obj.dofSpace.getDofIndices(obj.dofTypeIndices(4), elnodes);

                order = [1 3 9 7];
                X(el,:) = physics.mesh.Nodes(elnodes(order),1);
                Y(el,:) = physics.mesh.Nodes(elnodes(order),2);
                H(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,1));
				OH(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,2));%
				FE(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,5));
				FEOH(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,6));
				E(el,:) = physics.StateVec{obj.C_step}(Edofs(order));

				PHI(el,:) = physics.StateVec{obj.phi_step}(dofsPHI(order));
				CL(el,:) = physics.StateVec{obj.C_step}(dofsCL(order));
				T(el,:) = physics.StateVec{obj.C_step}(dofsT(order));
			end
            %patch(X',Y',Z','EdgeColor','None','FaceColor','interp');
			subplot(3,2,6)
			patch(X',Y',FE',FE','FaceColor','interp','EdgeColor','interp');
			title("FeOH")
			colorbar

			subplot(3,2,3)
			patch(X',Y',E',E','FaceColor','interp','EdgeColor','interp');
			title("E")
			colorbar

			subplot(3,2,1)
			H(H<=0)=NaN;
			H(~isnan(H)) = -log10(H(~isnan(H))/1000); %
			patch(X',Y',H',H','FaceColor','interp','EdgeColor','interp');
			title("pH")
			colorbar

			subplot(3,2,2)
			OH(OH<=0)=NaN;
			OH(~isnan(OH)) = -log10(OH(~isnan(OH))/1000); %
			patch(X',Y',OH',OH','FaceColor','interp','EdgeColor','interp');
			title("pOH")
			colorbar

			subplot(3,2,4)
			patch(X',Y',T',T','FaceColor','interp','EdgeColor','interp');
			title("\theta")
			colorbar

			subplot(3,2,5)
			patch(X',Y',PHI',PHI','FaceColor','interp','EdgeColor','interp');
			title("\varphi")
			colorbar
			
		end

		function plotFieldspH(obj, physics)
			for el=1:length(obj.ActiveElems)
                elnodes =physics.mesh.Elementgroups{obj.myGroupIndex}.Elems(obj.ActiveElems(el),:);
				Edofs = physics.dofSpace.getDofIndices(obj.fracTypeIndices(1), elnodes);
				for s=1:obj.n_species
					Cdofs(:,s) = physics.dofSpace.getDofIndices(obj.fracTypeIndices(2+s), elnodes);
				end
				dofsT = obj.dofSpace.getDofIndices(obj.fracTypeIndices(2), elnodes);
				dofsPHI= obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), elnodes);
				dofsCL= obj.dofSpace.getDofIndices(obj.dofTypeIndices(4), elnodes);

                order = [1 3 9 7];
                X(el,:) = physics.mesh.Nodes(elnodes(order),1);
                Y(el,:) = physics.mesh.Nodes(elnodes(order),2);
                H(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,1));
				OH(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,2));
				H(H<=0) = nan;
				pH(el,:) = -log10(H(el,:)/1000);
				FE(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,5));
				FEOH(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,6));
				E(el,:) = physics.StateVec{obj.C_step}(Edofs(order));

				PHI(el,:) = physics.StateVec{obj.phi_step}(dofsPHI(order));
				CL(el,:) = physics.StateVec{obj.C_step}(dofsCL(order));
				T(el,:) = physics.StateVec{obj.C_step}(dofsT(order));
			end
			pH(PHI<0.1) = nan;
			patch(X',Y',pH',pH','FaceColor','interp','EdgeColor','interp');
		end

		function plotFieldsFe(obj, physics)
			for el=1:length(obj.ActiveElems)
                elnodes =physics.mesh.Elementgroups{obj.myGroupIndex}.Elems(obj.ActiveElems(el),:);
				Edofs = physics.dofSpace.getDofIndices(obj.fracTypeIndices(1), elnodes);
				for s=1:obj.n_species
					Cdofs(:,s) = physics.dofSpace.getDofIndices(obj.fracTypeIndices(2+s), elnodes);
				end
				dofsT = obj.dofSpace.getDofIndices(obj.fracTypeIndices(2), elnodes);
				dofsPHI= obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), elnodes);
				dofsCL= obj.dofSpace.getDofIndices(obj.dofTypeIndices(4), elnodes);

                order = [1 3 9 7];
                X(el,:) = physics.mesh.Nodes(elnodes(order),1);
                Y(el,:) = physics.mesh.Nodes(elnodes(order),2);
                H(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,1));
				OH(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,2));
				H(H<=0) = nan;
				pH(el,:) = -log10(H(el,:)/1000);
				FE(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,5));
				FEOH(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,6));
				E(el,:) = physics.StateVec{obj.C_step}(Edofs(order));

				PHI(el,:) = physics.StateVec{obj.phi_step}(dofsPHI(order));
				CL(el,:) = physics.StateVec{obj.C_step}(dofsCL(order));
				T(el,:) = physics.StateVec{obj.C_step}(dofsT(order));
			end
			FE(PHI<0.1) = nan;
			patch(X',Y',FE',FE','FaceColor','interp','EdgeColor','interp');
		end
        
		function plotFieldsSurf(obj, physics)
			for el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                elnodes =physics.mesh.Elementgroups{obj.myGroupIndex}.Elems(el,:);
				Edofs = physics.dofSpace.getDofIndices(obj.dofTypeIndices(1), elnodes);
				for s=1:obj.n_species
					Cdofs(:,s) = physics.dofSpace.getDofIndices(obj.dofTypeIndices(1+s), elnodes);
				end

                order = [1 2 3]; % [1 4 2 5 3 6]
                X(el,:) = physics.mesh.Nodes(elnodes(order),1);
                Y(el,:) = physics.mesh.Nodes(elnodes(order),2);
                H(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,1));
				OH(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,2));%
				FE(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,5));
				FEOH(el,:) = physics.StateVec{obj.C_step}(Cdofs(order,6));
				E(el,:) = physics.StateVec{obj.C_step}(Edofs(order));
			end
            %patch(X',Y',Z','EdgeColor','None','FaceColor','interp');
			subplot(2,2,1)
			patch(X',Y',H',H','FaceColor','interp','EdgeColor','interp');
			title("C_{H^+}")

			subplot(2,2,2)
			patch(X',Y',OH',OH','FaceColor','interp','EdgeColor','interp');
			title("C_{OH^-}")
				
			subplot(2,2,3)
			patch(X',Y',FE',FE','FaceColor','interp','EdgeColor','interp');
			title("FeOH^+")

			subplot(2,2,4)
			patch(X',Y',E',E','FaceColor','interp','EdgeColor','interp');
			title("E")
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
           
        end
        
    end
end

