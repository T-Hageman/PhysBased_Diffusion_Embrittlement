classdef HydrogenDiffusionUnDamaged < BaseModel
    %HYDROGENDIFFUSION Physics model implementing the conservation of
	%species within the metal interstitial lattice sites, assuming the material is UNDAMAGED
	% , including diffusivity due to concentration gradients and hydrostatic stress
	%gradients (NOT taking into account the current phase-field state). 
	%Input properties required:
	% physics_in{2}.type = "HydrogenDiffusionUnDamaged";
	% physics_in{2}.Egroup = "Internal";
    % physics_in{2}.DL = 1e-9;				%Diffusivity [m/s]
	% physics_in{2}.NL = physics_in{1}.NL;	%Concentration of interstitial lattice sites [mol/m^3]	
	% physics_in{2}.gb = physics_in{1}.gb;	%Grain boundary binding energy [J/mol]
	% physics_in{2}.NT = 1e2;					%concentration of trapping sites
    
    properties
        mesh			%Pointer to mesh object
        myName			%Name of this model
        myGroup			%String indicating the element group this model operates on
        myGroupIndex	%Index of the element group involved in this model
        dofSpace		%pointer to degree of freedom space
        dofTypeIndices	%Indices of degrees of freedom required by this model
        
        DL			%Lattice diffusivity [m/s]
		NL			%Interstitial lattice site concentration [mol/m^3]
        poisson		%Poisson ratio of the metal [-]
        young		%Youngs modulus of the metal [Pa]

		gb			%binding energy of trapping sites
		NT			%concentration of trapping sites [mol/m^3]

		CL_int		%Current total hydrogen contents of the domain
		CL_max		%Current maximum hydrogen concentration within the domain

		R_const = 8.31446261815324;	%gas constant
		T_const = 293.15;			%reference temperature
		VH_const = 2e-6;			%hydrogen volume

		dx_Step		%Step in which displacements are resolved
		C_Step		%step in which chemical sub-problem is resolved
    end
    
    methods
        function obj = HydrogenDiffusionUnDamaged(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "HydrogenDiffusion";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            [obj.dofTypeIndices, stp] = obj.dofSpace.getDofType({"dx","dy","CL"});
			obj.dx_Step = stp(1);
			obj.C_Step  = stp(3);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            %% get parameters
            obj.DL = inputs.DL;
			obj.NL = inputs.NL;
			obj.gb = inputs.gb;
			obj.NT = inputs.NT;
            
			obj.poisson = physics.models{1}.poisson;
            obj.young = physics.models{1}.young;
        end
        
        function getKf(obj, physics, stp)
			if (stp == obj.C_Step)
            	fprintf("        HydrogenDiffusion get Matrix:")
            	t = tic;
            	
            	dt = physics.dt;
				Svec = physics.StateVec;
            	SvecOld = physics.StateVec_Old;
	
            	dofmatX = [];
            	dofmatY = [];
            	kmat = [];
				fvec = [];
            	dofvec = [];
	
				CL_sum = 0;
				CL_max2 = 0;
				parfor n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                	Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                	[N, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
					G2 = obj.mesh.getG2(obj.myGroupIndex, n_el);
	
                	dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                	dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                	dofsCL= obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
                	dofsXY = [dofsX; dofsY];
	
                	X = Svec{obj.dx_Step}(dofsX);
                	Y = Svec{obj.dx_Step}(dofsY);
                	XY = [X;Y];
                	CL = Svec{obj.C_Step}(dofsCL);
                	CLOld = SvecOld{obj.C_Step}(dofsCL);
	
                	q_el = zeros(length(dofsCL), 1);
	
                	K_cu = zeros(length(dofsCL), length(dofsXY));
                	K_cc = zeros(length(dofsCL));
					for ip=1:length(w)
                    	%% capacity term
                    	q_el = q_el + w(ip) * N(ip,:)'*N(ip,:)*(CL-CLOld)/dt;
                    	K_cc = K_cc + w(ip) * N(ip,:)'*N(ip,:)/dt;

						%% capacity term (traps)
						dNTdC = obj.NT/obj.NL*exp(-obj.gb/obj.R_const/obj.T_const)/((max(0,N(ip,:)*CL)/obj.NL+exp(-obj.gb/obj.R_const/obj.T_const))^2);
						dNTdC_dC = -2*obj.NT/obj.NL^2*exp(-obj.gb/obj.R_const/obj.T_const)/((max(0,N(ip,:)*CL)/obj.NL+exp(-obj.gb/obj.R_const/obj.T_const))^3);
						q_el = q_el + w(ip) * N(ip,:)'*dNTdC*N(ip,:)*(CL-CLOld)/dt;
						K_cc = K_cc + w(ip) * N(ip,:)'*dNTdC*N(ip,:)/dt + 0*N(ip,:)'*dNTdC_dC*N(ip,:)*(N(ip,:)*(CL-CLOld))/dt;
	
                    	%% hydraulic stress driven  
						pfx = obj.young/(3*(1-2*obj.poisson));
                    	Bstar = obj.getBstar(G2(ip,:,:));
						BLine = [squeeze(G(ip,:,1)) squeeze(G(ip,:,2))];
                    	dsh = pfx*Bstar*XY;
	
 						q_el = q_el - w(ip)*obj.DL*obj.VH_const/obj.R_const/obj.T_const * squeeze(G(ip,:,:)) * dsh *max(0,(N(ip,:)*CL));
 						K_cc = K_cc - w(ip)*obj.DL*obj.VH_const/obj.R_const/obj.T_const * (squeeze(G(ip,:,:)) * dsh) *N(ip,:);
						%K_cu = K_cu - w(ip)*obj.DL*obj.VH_const/obj.R_const/obj.T_const * squeeze(G(ip,:,:)) * pfx*Bstar*max(0,(N(ip,:)*CL));
	
                    	%% diffusion driven
                    	q_el = q_el + w(ip)*obj.DL*(1/max(1e-20,(1-max(0,N(ip,:)*CL)/obj.NL)))*squeeze(G(ip,:,:))*squeeze(G(ip,:,:))'*CL;
                    	K_cc = K_cc + w(ip)*obj.DL*(1/max(1e-20,(1-max(0,N(ip,:)*CL)/obj.NL)))*squeeze(G(ip,:,:))*squeeze(G(ip,:,:))' ....
						        	+ w(ip)*obj.DL*(1/max(1e-20,(1-max(0,N(ip,:)*CL)/obj.NL))^2)*squeeze(G(ip,:,:))*((squeeze(G(ip,:,:))'*CL)*N(ip,:))/obj.NL;  %THIS TERM WAS *0
	
						CL_sum = CL_sum + w(ip)*N(ip,:)*CL;
						CL_max2 = max(CL_max2, N(ip,:)*CL);
					end
	
%                 	[dofmatxloc,dofmatyloc] = ndgrid(dofsCL,dofsXY);
%                 	dofmatX = [dofmatX; dofmatxloc(:)];
%                 	dofmatY = [dofmatY; dofmatyloc(:)];
%                 	kmat = [kmat; K_cu(:)];
	
                	[dofmatxloc,dofmatyloc] = ndgrid(dofsCL,dofsCL);
                	dofmatX = [dofmatX; dofmatxloc(:)];
                	dofmatY = [dofmatY; dofmatyloc(:)];
                	kmat = [kmat; K_cc(:)];        
	
					fvec = [fvec; q_el];
                	dofvec = [dofvec; dofsCL];
				end 
	
				physics.fint{obj.C_Step} = physics.fint{obj.C_Step} + sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint{obj.C_Step}), 1);
            	physics.K{obj.C_Step} = physics.K{obj.C_Step} + sparse(dofmatX, dofmatY, kmat, length(physics.fint{obj.C_Step}),length(physics.fint{obj.C_Step}));
					
				obj.CL_int = CL_sum;
				obj.CL_max = CL_max2;
            	
            	tElapsed = toc(t);
            	fprintf("            (Assemble time:"+string(tElapsed)+")\n");
			end
        end
        
		function B = getBstar(~, grad2s)
            cp_count = size(grad2s, 2);
            B = zeros(2, cp_count*2);
            for ii = 1:cp_count %using plane strain e_zz = 0
				%dx
				B(1, ii) = grad2s(1,ii, 1);
				B(2, ii) = grad2s(1,ii, 3);

				%dy
				B(1, ii + cp_count) = grad2s(1,ii, 3);
				B(2, ii + cp_count) = grad2s(1,ii, 2);
            end
        end

        function [hasInfo, provided] = Provide_Info(obj, physics, var, elems, loc)
           hasInfo = false;
           provided = [];
           
        end
        
    end
end

