classdef ElectrolyteInterface < BaseModel
    %ELECTROLYTEINTERFACE Physics model implementing surface reactions at
	%the metal-electrolyte interface, and implementing the mass balance for
	%the surface adsorbed hydrogen. The surface reactions resolved are the
	%individual steps for the hydrogen evolution reaction, and the
	%corrosion reaction (assuming neglible surface dissolution, and thus 
	%not geometryically represented). The input parameters required for
	%this model are:
	%	physics_in{9}.type = "ElectrolyteInterface";
    %	physics_in{9}.Egroup = "Interface";
	%	physics_in{9}.NAds = 1e-3; % Concentration of surface sites [mol/m^2]
	%	physics_in{9}.k = [	1e-4,	1e-10,	0.5,	0;   %reaction constants, [k, k', alpha, E_eq] for	Acidic Volmer
	%						1e-10,	0,		0.3,	0;	 %											Acidic Heyrovsky
	%						1e-6,	0,		0,		0;	 %											Tafel
	%						1e3,	7e7,	0,		0;	 %											Absorbtion
	%						1e-8,	1e-13,	0.5,	0;	 %											Basic Volmer
	%						1e-10,	1e-14,	0.3,	0;   %											Basic Heyrovsky
	%						3e-5/(2*96485.3329),3e-5/(2*96485.3329), 0.5, -0.4];   %				Corrosion
	%	physics_in{9}.NL = physics_in{2}.NL; % Concentration of interstitial lattice sites [mol/m^3]
	%	physics_in{9}.Em = 0.0; % Metal Potential [V_SHE]
	%	physics_in{9}.Lumped = [1 1 1 1 1 1 1]; %Flags to enable lumped integration on a per-reaction basis
    
    properties
        mesh			%Pointer to mesh object
        myName			%Name of this model (for identification purposes only)
        myGroup			%String indicating which element group this model operates on
        myGroupIndex	%Index of element group involved with this model
        dofSpace		%Pointer to degree of freedom management object
        dofTypeIndices	%Indices of degrees of freedom

		k				%Reaction rate matrix, [k, k', alpha, E_eq], 1 row per step in order: acidic volmer, acidic heyrovsky, tafel, absorbtion, basic volmer, basic heyrovsky, corrosion
		Em				%Applied metal potential [V_SHE]
		NL				%Concentration of interstitial lattice sites [mol/m^3]
		NAds			%Concentration of surface adsorption sites [mol/m^2]
		Lumped			%Flags to enable lumped inegration on a per-reaction basis

		C_Step

		F_const = 96485.3329;		%Faraday constant
		R_const = 8.31446261815324;	%Gas constant
		T_const = 293.15;			%Reference temperature

		n_species
    end
    
    methods
        function obj = ElectrolyteInterface(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "ElectrolyteInterface";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs	                       %1       2      3    4     5    6    7     8      9
            [obj.dofTypeIndices, stp] = obj.dofSpace.getDofType({"Epot", "theta","CL","H", "OH","Na", "Cl","Fe","FeOH"});
			obj.C_Step = stp(1);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            %% get parameters
			obj.Lumped = inputs.Lumped;
			obj.k = inputs.k;

			obj.NAds = inputs.NAds;
			obj.Em = inputs.Em;
			obj.NL = inputs.NL;

			obj.n_species = 6;
        end
        
        function getKf(obj, physics, stp)
			if (obj.C_Step == stp)
            	fprintf("        ElectrolyteInterface get Matrix:")
            	t = tic;
            	
            	dt = physics.dt;
	
				%stiffness matrix
            	dofmatX = [];
            	dofmatY = [];
            	kmat = [];
            	fvec = [];
            	dofvec = [];
	
				Svec = physics.StateVec;
            	SvecOld = physics.StateVec_Old;

            	parfor n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                	Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                	[N, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
	
					dofsE = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
					dofsT = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
					dofsC = zeros(length(dofsE), obj.n_species);
					for s=1:obj.n_species
						dofsC(:,s) = obj.dofSpace.getDofIndices(obj.dofTypeIndices(s+3), Elem_Nodes);
					end
	
					dofsCL= obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
	
					C = Svec{obj.C_Step}(dofsC);
					E = Svec{obj.C_Step}(dofsE);
					T = Svec{obj.C_Step}(dofsT);
					CL = Svec{obj.C_Step}(dofsCL);
	
					C_OLD = SvecOld{obj.C_Step}(dofsC);
					E_OLD = SvecOld{obj.C_Step}(dofsE);
					T_OLD = SvecOld{obj.C_Step}(dofsT);
					CL_Old = SvecOld{obj.C_Step}(dofsCL);
	
					q_C    = zeros(length(dofsE), obj.n_species);
					dqC_dC = zeros(length(dofsE), length(dofsE), obj.n_species, obj.n_species);
					dqC_dE = zeros(length(dofsE), length(dofsE), obj.n_species);
					dqC_dT = zeros(length(dofsE), length(dofsE), obj.n_species);
	
					q_T = zeros(length(dofsE), 1);
					dqT_dT = zeros(length(dofsE), length(dofsE));
					dqT_dE = zeros(length(dofsE), length(dofsE));
					dqT_dC = zeros(length(dofsE), length(dofsE), obj.n_species);
					dqT_dCL = zeros(length(dofsE), length(dofsCL));
	
					q_CL = zeros(length(dofsCL), 1);
					dqCL_dCL = zeros(length(dofsCL), length(dofsCL));
					dqCL_dT = zeros(length(dofsCL), length(dofsE));
	
					C_Lumped = zeros(length(dofsE), 1);
					for ip=1:length(w)
	
						% surface capacity
						q_T = q_T + w(ip)*obj.NAds/dt*N(ip,:)'*N(ip,:)*(T-T_OLD);
						dqT_dT = dqT_dT+ w(ip)*obj.NAds/dt*N(ip,:)'*N(ip,:);
	
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
							q_C(:,1) = q_C(:,1) - w(ip)*N(ip,:)'*(react(r,1)-react(r,2))*products(r,1) *(1-obj.Lumped(r));
							q_C(:,2) = q_C(:,2) - w(ip)*N(ip,:)'*(react(r,1)-react(r,2))*products(r,2)*(1-obj.Lumped(r));
							q_C(:,5) = q_C(:,5) - w(ip)*N(ip,:)'*(react(r,1)-react(r,2))*products(r,3)*(1-obj.Lumped(r));
							q_T      = q_T      - w(ip)*N(ip,:)'*(react(r,1)-react(r,2))*products(r,4) *(1-obj.Lumped(r));
							q_CL     = q_CL     - w(ip)*N(ip,:)'*(react(r,1)-react(r,2))*products(r,5)*(1-obj.Lumped(r));
	
							for n=1:obj.n_species
								dqC_dC(:,:,1,n) = dqC_dC(:,:,1,n)-w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,1)*(1-obj.Lumped(r));
								dqC_dC(:,:,2,n) = dqC_dC(:,:,2,n)-w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,2)*(1-obj.Lumped(r));
								dqC_dC(:,:,5,n) = dqC_dC(:,:,5,n)-w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,3)*(1-obj.Lumped(r));
								dqT_dC(:,:,n)   = dqT_dC(:,:,n)  -w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,4)*(1-obj.Lumped(r));
							end
							dqC_dE(:,:,1)   = dqC_dE(:,:,1) - w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,1)-dreact(r,2,1))*products(r,1)*(1-obj.Lumped(r));
							dqC_dE(:,:,2)   = dqC_dE(:,:,2) - w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,1)-dreact(r,2,1))*products(r,2)*(1-obj.Lumped(r));
							dqC_dE(:,:,5)   = dqC_dE(:,:,5) - w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,1)-dreact(r,2,1))*products(r,3)*(1-obj.Lumped(r));
							dqT_dE          = dqT_dE        - w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,1)-dreact(r,2,1))*products(r,4)*(1-obj.Lumped(r));
	
							dqC_dT(:,:,1)   = dqC_dT(:,:,1) - w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,2)-dreact(r,2,2))*products(r,1)*(1-obj.Lumped(r));
							dqC_dT(:,:,2)   = dqC_dT(:,:,2) - w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,2)-dreact(r,2,2))*products(r,2)*(1-obj.Lumped(r));
							dqC_dT(:,:,5)   = dqC_dT(:,:,5) - w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,2)-dreact(r,2,2))*products(r,3)*(1-obj.Lumped(r));
							dqT_dT          = dqT_dT        - w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,2)-dreact(r,2,2))*products(r,4)*(1-obj.Lumped(r));
							dqCL_dT         = dqCL_dT       - w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,2)-dreact(r,2,2))*products(r,5)*(1-obj.Lumped(r));
	
							dqT_dCL         = dqT_dCL       - w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,3)-dreact(r,2,3))*products(r,4)*(1-obj.Lumped(r));
							dqCL_dCL        = dqCL_dCL      - w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,3)-dreact(r,2,3))*products(r,5)*(1-obj.Lumped(r));
						end
	
						% lumped integration weight
						C_Lumped = C_Lumped + w(ip)*N(ip,:)';
					end
	
	
					for i=1:length(dofsT) %% lumped integrations
						[react, dreact, products] = obj.reactions(C(i,1), C(i,2),C(i,5), T(i), E(i), CL(i));
	
						for r=1:7
							q_C(i,1) = q_C(i,1) - C_Lumped(i)*(react(r,1)-react(r,2))*products(r,1)*obj.Lumped(r);
							q_C(i,2) = q_C(i,2) - C_Lumped(i)*(react(r,1)-react(r,2))*products(r,2)*obj.Lumped(r);
							q_C(i,5) = q_C(i,5) - C_Lumped(i)*(react(r,1)-react(r,2))*products(r,3)*obj.Lumped(r);
							q_T(i)   = q_T(i)   - C_Lumped(i)*(react(r,1)-react(r,2))*products(r,4)*obj.Lumped(r);
							q_CL(i)  = q_CL(i)  - C_Lumped(i)*(react(r,1)-react(r,2))*products(r,5)*obj.Lumped(r);
	
							for n=1:obj.n_species
								dqC_dC(i,i,1,n) = dqC_dC(i,i,1,n)-C_Lumped(i)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,1)*obj.Lumped(r);
								dqC_dC(i,i,2,n) = dqC_dC(i,i,2,n)-C_Lumped(i)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,2)*obj.Lumped(r);
								dqC_dC(i,i,5,n) = dqC_dC(i,i,5,n)-C_Lumped(i)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,3)*obj.Lumped(r);
								dqT_dC(i,i,n)   = dqT_dC(i,i,n)  -C_Lumped(i)*(dreact(r,1,3+n)-dreact(r,2,3+n))*products(r,4)*obj.Lumped(r);
							end
							dqC_dE(i,i,1)   = dqC_dE(i,i,1) - C_Lumped(i)*(dreact(r,1,1)-dreact(r,2,1))*products(r,1)*obj.Lumped(r);
							dqC_dE(i,i,2)   = dqC_dE(i,i,2) - C_Lumped(i)*(dreact(r,1,1)-dreact(r,2,1))*products(r,2)*obj.Lumped(r);
							dqC_dE(i,i,5)   = dqC_dE(i,i,5) - C_Lumped(i)*(dreact(r,1,1)-dreact(r,2,1))*products(r,3)*obj.Lumped(r);
							dqT_dE(i,i)     = dqT_dE(i,i)   - C_Lumped(i)*(dreact(r,1,1)-dreact(r,2,1))*products(r,4)*obj.Lumped(r);
	
							dqC_dT(i,i,1)   = dqC_dT(i,i,1) - C_Lumped(i)*(dreact(r,1,2)-dreact(r,2,2))*products(r,1)*obj.Lumped(r);
							dqC_dT(i,i,2)   = dqC_dT(i,i,2) - C_Lumped(i)*(dreact(r,1,2)-dreact(r,2,2))*products(r,2)*obj.Lumped(r);
							dqC_dT(i,i,5)   = dqC_dT(i,i,5) - C_Lumped(i)*(dreact(r,1,2)-dreact(r,2,2))*products(r,3)*obj.Lumped(r);
							dqT_dT(i,i)     = dqT_dT(i,i)   - C_Lumped(i)*(dreact(r,1,2)-dreact(r,2,2))*products(r,4)*obj.Lumped(r);
							dqCL_dT(i,i)    = dqCL_dT(i,i)  - C_Lumped(i)*(dreact(r,1,2)-dreact(r,2,2))*products(r,5)*obj.Lumped(r);
	
							dqT_dCL(i,i)    = dqT_dCL(i,i)  - C_Lumped(i)*(dreact(r,1,3)-dreact(r,2,3))*products(r,4)*obj.Lumped(r);
							dqCL_dCL(i,i)   = dqCL_dCL(i,i) - C_Lumped(i)*(dreact(r,1,3)-dreact(r,2,3))*products(r,5)*obj.Lumped(r);
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
						tmp = dqC_dE(:,:,s1);
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
	
					fvec = [fvec; q_T];
                	dofvec = [dofvec; dofsT];
	
                	fvec = [fvec; q_CL];
                	dofvec = [dofvec; dofsCL];
            	end 
	
            	physics.fint{obj.C_Step} = physics.fint{obj.C_Step} + sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint{obj.C_Step}), 1);
            	physics.K{obj.C_Step} = physics.K{obj.C_Step} + sparse(dofmatX, dofmatY, kmat, length(physics.fint{obj.C_Step}),length(physics.fint{obj.C_Step}));
	
				%fprintf("\n"+string(fraccheck)+"\n");
	
				%figure(753)
				%clf
            	%obj.plotFields(physics)
            	tElapsed = toc(t);
            	fprintf("            (Assemble time:"+string(tElapsed)+")\n");
			end
		end

		function plotFields(obj, physics)
			for el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                elnodes =physics.mesh.Elementgroups{obj.myGroupIndex}.Elems(el,:);
				Edofs = physics.dofSpace.getDofIndices(obj.dofTypeIndices(1), elnodes);
				Tdofs = physics.dofSpace.getDofIndices(obj.dofTypeIndices(2), elnodes);
				for s=1:obj.n_species
					Cdofs(:,s) = physics.dofSpace.getDofIndices(obj.dofTypeIndices(3+s), elnodes);
				end

                order = [1 3 2];
                X(el,:) = physics.mesh.Nodes(elnodes(order),1);
                Y(el,:) = physics.mesh.Nodes(elnodes(order),2);
                H(el,:) = -log10(physics.StateVec{obj.C_Step}(Cdofs(order,1))/1000);
				OH(el,:) = -log10(physics.StateVec{obj.C_Step}(Cdofs(order,2))/1000);
				FE(el,:) = physics.StateVec{obj.C_Step}(Cdofs(order,5));
				FEOH(el,:) = physics.StateVec{obj.C_Step}(Cdofs(order,6));
				E(el,:) = physics.StateVec{obj.C_Step}(Edofs(order));
				T(el,:) = physics.StateVec{obj.C_Step}(Tdofs(order));
			end
            %patch(X',Y',Z','EdgeColor','None','FaceColor','interp');
			subplot(2,3,1)
			fill3(X',Y',H',H','FaceColor','interp');
			title("pH")

			subplot(2,3,2)
			fill3(X',Y',OH',OH','FaceColor','interp');
			title("pOH")
				
			subplot(2,3,3)
			fill3(X',Y',FEOH',FEOH','FaceColor','interp');
			title("FeOH")

			subplot(2,3,4)
			fill3(X',Y',E',E','FaceColor','interp');
			title("E")

			subplot(2,3,5)
			fill3(X',Y',T',T','FaceColor','interp');
			title("T")

			subplot(2,3,6)
			fill3(X',Y',FE',FE','FaceColor','interp');
			title("Fe")


		end

		function plotReactions(obj, physics)
			for el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                elnodes =physics.mesh.Elementgroups{obj.myGroupIndex}.Elems(el,:);
				Edofs = physics.dofSpace.getDofIndices(obj.dofTypeIndices(1), elnodes);
				Tdofs = physics.dofSpace.getDofIndices(obj.dofTypeIndices(2), elnodes);
				dofsCL= physics.dofSpace.getDofIndices(obj.dofTypeIndices(3), elnodes);
				for s=1:obj.n_species
					Cdofs(:,s) = physics.dofSpace.getDofIndices(obj.dofTypeIndices(3+s), elnodes);
				end

                order = [1 2 3];
                X(el,:) = [physics.mesh.Nodes(elnodes(order),1);NaN];
                Y(el,:) = [physics.mesh.Nodes(elnodes(order),2);NaN];				

				for i=1:length(order)
					CH = physics.StateVec{obj.C_Step}(Cdofs(order(i),1));
					COH= physics.StateVec{obj.C_Step}(Cdofs(order(i),2));
					CFE = physics.StateVec{obj.C_Step}(Cdofs(order(i),5));
					theta = physics.StateVec{obj.C_Step}(Tdofs(order(i)));
					phil = physics.StateVec{obj.C_Step}(Edofs(order(i)));
					CLat = physics.StateVec{obj.C_Step}(dofsCL(order(i)));

					[react, ~, ~] = obj.reactions(CH, COH, CFE, theta, phil, CLat);
					r(el,i,:) = react(:,1)-react(:,2);
				end
				r(el,4,:) = NaN;
			end
            %patch(X',Y',Z','EdgeColor','None','FaceColor','interp');
			subplot(3,3,1)
			fill3(X',Y',r(:,:,1)',r(:,:,1)','FaceColor','interp');
			title("\nu_1")

			subplot(3,3,2)
			fill3(X',Y',r(:,:,2)',r(:,:,2)','FaceColor','interp');
			title("\nu_2")
				
			subplot(3,3,3)
			fill3(X',Y',r(:,:,3)',r(:,:,3)','FaceColor','interp');
			title("\nu_3")

			subplot(3,3,4)
			fill3(X',Y',r(:,:,4)',r(:,:,4)','FaceColor','interp');
			title("\nu_4")

			subplot(3,3,5)
			fill3(X',Y',r(:,:,5)',r(:,:,5)','FaceColor','interp');
			title("\nu_5")

			subplot(3,3,6)
			fill3(X',Y',r(:,:,6)',r(:,:,6)','FaceColor','interp');
			title("\nu_6")

			subplot(3,3,7)
			fill3(X',Y',r(:,:,7)',r(:,:,4)','FaceColor','interp');
			title("\nu_7")


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

					dCFE = 1;
					if (CFE>Cmax)
						CFE = Cmax;
						dCFE = 0;
					end

					% reaction 1 
					react(1,1)    = obj.k(1,1)*CH *(1-theta)*exp(-obj.k(1,3)*(obj.Em-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(1,1,2) = obj.k(1,1)*CH *-1              *exp(-obj.k(1,3)*(obj.Em-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(1,1,4) = obj.k(1,1)*dCH*(1-theta)*exp(-obj.k(1,3)*(obj.Em-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(1,1,1) = obj.k(1,1)*CH *(1-theta)*exp(-obj.k(1,3)*(obj.Em-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.k(1,3)*(-1)*obj.F_const/obj.R_const/obj.T_const);

					react(1,2)    = obj.k(1,2)*theta *exp((1-obj.k(1,3))*(obj.Em-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(1,2,2) = obj.k(1,2)              *exp((1-obj.k(1,3))*(obj.Em-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(1,2,1) = obj.k(1,2)*theta *exp((1-obj.k(1,3))*(obj.Em-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const)*((1-obj.k(1,3))*(-1)*obj.F_const/obj.R_const/obj.T_const);

					% reaction 2
					react(2,1)    = obj.k(2,1)*CH *theta*exp(-obj.k(2,3)*(obj.Em-phil-obj.k(2,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(2,1,4) = obj.k(2,1)*dCH*theta*exp(-obj.k(2,3)*(obj.Em-phil-obj.k(2,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(2,1,2) = obj.k(2,1)*CH              *exp(-obj.k(2,3)*(obj.Em-phil-obj.k(2,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(2,1,1) = obj.k(2,1)*CH *theta*exp(-obj.k(2,3)*(obj.Em-phil-obj.k(2,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.k(2,3)*(-1)*obj.F_const/obj.R_const/obj.T_const);

					%react(2,2) = 0;

					%reaction 3
					react(3,1)    = obj.k(3,1)*abs(theta)*theta;
					dreact(3,1,2) = 2*obj.k(3,1)*abs(theta);

					%react(3,2) = 0;

					%reaction 4
					%if theta>=0
						react(4,1)    = obj.k(4,1)*(obj.NL-max(0,CLat))*theta;
						dreact(4,1,2) = obj.k(4,1)*(obj.NL-max(0,CLat));
						dreact(4,1,3) = obj.k(4,1)*(-1)*theta;
					%end
					react(4,2)    =  obj.k(4,2)*CLat*(1-theta);
					dreact(4,2,2) =  obj.k(4,2)*CLat*(-1);
					dreact(4,2,3) =  obj.k(4,2)*(1-theta);


					%reaction 5
					react(5,1)    = obj.k(5,1)*(1-theta)*exp(-obj.k(5,3)*(obj.Em-phil-obj.k(5,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(5,1,2) = obj.k(5,1)*-1       *exp(-obj.k(5,3)*(obj.Em-phil-obj.k(5,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(5,1,1) = obj.k(5,1)*(1-theta)*exp(-obj.k(5,3)*(obj.Em-phil-obj.k(5,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.k(5,3)*(-1)*obj.F_const/obj.R_const/obj.T_const);

					react(5,2)    = obj.k(5,2)*COH *theta*exp((1-obj.k(5,3))*(obj.Em-phil-obj.k(5,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(5,2,2) = obj.k(5,2)*COH       *exp((1-obj.k(5,3))*(obj.Em-phil-obj.k(5,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(5,2,1) = obj.k(5,2)*COH *theta*exp((1-obj.k(5,3))*(obj.Em-phil-obj.k(5,4))*obj.F_const/obj.R_const/obj.T_const)*((1-obj.k(5,3))*(-1)*obj.F_const/obj.R_const/obj.T_const);
					dreact(5,2,5) = obj.k(5,2)*dCOH*theta*exp((1-obj.k(5,3))*(obj.Em-phil-obj.k(5,4))*obj.F_const/obj.R_const/obj.T_const);

					%reaction 6
					react(6,1)    = obj.k(6,1)*theta*exp(-obj.k(6,3)*(obj.Em-phil-obj.k(6,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(6,1,2) = obj.k(6,1)      *exp(-obj.k(6,3)*(obj.Em-phil-obj.k(6,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(6,1,1) = obj.k(6,1)*theta*exp(-obj.k(6,3)*(obj.Em-phil-obj.k(6,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.k(6,3)*(-1)*obj.F_const/obj.R_const/obj.T_const);

					%react(6,2) = 0

					%corrosion
% 					react(7,1)    = obj.k(7,1)*(1-CFE/Cmax)*exp((1-obj.k(7,3))*(obj.Em-phil-obj.k(7,4))*obj.F_const/obj.R_const/obj.T_const);
% 					dreact(7,1,1) = obj.k(7,1)*(1-CFE/Cmax)*exp((1-obj.k(7,3))*(obj.Em-phil-obj.k(7,4))*obj.F_const/obj.R_const/obj.T_const)*((1-obj.k(7,3))*(-1)*obj.F_const/obj.R_const/obj.T_const);
% 					dreact(7,1,8) = obj.k(7,1)*(-dCFE/Cmax)*exp((1-obj.k(7,3))*(obj.Em-phil-obj.k(7,4))*obj.F_const/obj.R_const/obj.T_const);


					react(7,1)    = obj.k(7,1)*CFE*exp(-obj.k(7,3)*(obj.Em-phil-obj.k(7,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(7,1,1) = obj.k(7,1)*CFE*exp(-obj.k(7,3)*(obj.Em-phil-obj.k(7,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.k(7,3)*(-1)*obj.F_const/obj.R_const/obj.T_const);
					dreact(7,1,8) = obj.k(7,1)*exp(-obj.k(7,3)*(obj.Em-phil-obj.k(7,4))*obj.F_const/obj.R_const/obj.T_const);
			
					react(7,2)    = obj.k(7,2)*exp((1-obj.k(7,3))*(obj.Em-phil-obj.k(7,4))*obj.F_const/obj.R_const/obj.T_const);
					dreact(7,2,1) = obj.k(7,2)*exp((1-obj.k(7,3))*(obj.Em-phil-obj.k(7,4))*obj.F_const/obj.R_const/obj.T_const)*((1-obj.k(7,3))*(-1)*obj.F_const/obj.R_const/obj.T_const);
			
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

