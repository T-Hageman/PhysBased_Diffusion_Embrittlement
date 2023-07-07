classdef Physics < handle
    %PHYSICS Class that handles the definitions and assembly of stiffness
	%matrices, constraints, and statevectors
    
    properties
        mesh		%pointer to the mesh object
        models		%array containing all the included physical models
        dofSpace	%pointer to the degree of freedom object
        
        K			%Tangential stiffness matrix
        fint		%Internal force vector
        StateVec	%State vectors which is iteratively updated to obtain results at t+dt
        StateVec_Old	%Converged state vectors at time t
        time		%current time within the simulation
        dt			%Time increment size
        
        condofs		%constrained degrees of freedom
        convals		%values applied as constraints to the degrees of freedom
        convals_corr	%increment required to maintain constraint
        conMat		%matrix to transfer from unconstrained system to solely the constrained dofs
        unconMat	%matrix to transfer from unconstrained system to constrained
		nonz		%number of non-zero values in stiffness matrix
    end
    
    methods
        PlotNodal(obj, dofName, dispscale, plotloc)	%exterior defined, plots nodal quantities
        PlotIP(obj, varName, plotloc)				%exterior defined, plots integration point quantities
        
        function obj = Physics(mesh, inputs, dofs_in)
			%initialization

            obj.mesh = mesh;
            obj.dofSpace = DofSpace(obj.mesh, dofs_in);
            
            for i=1:length(inputs)
                f = str2func(inputs{i}.type);
                obj.models{i} = f(mesh, obj, inputs{i});
            end
            
            obj.dt = 0;
            
			for i=1:obj.dofSpace.NSteps
            	dofcount = obj.dofSpace.NDofs(i);
            	obj.StateVec{i} = zeros(dofcount, 1);
            	obj.StateVec_Old{i} = obj.StateVec{i};
            	
				obj.nonz(i) = 0;
            	obj.K{i} = [];
            	obj.fint{i} = zeros(dofcount,1);
			end

        end
        
		function OncePerStep(obj, stp)
			%procedures that should be performed once per step

            for m=1:length(obj.models)
                obj.models{m}.OncePerStep(obj, stp);
			end
		end


        function Assemble(obj, stp)
			%Assemble stiffness matrix and internal force vector for the
			%current step

            dofcount = obj.dofSpace.NDofs(stp);

			obj.condofs{stp} = [];
			obj.convals{stp} = [];
			
			if isempty(obj.K{stp})

			else
				obj.nonz(stp) = round(nnz(obj.K{stp}));
			end
            obj.K{stp} = spalloc(dofcount, dofcount, obj.nonz(stp));
            obj.fint{stp} = zeros(dofcount, 1);

            disp("    Assembling:")
            for m=1:length(obj.models)
                obj.models{m}.getKf(obj, stp);
            end
        end
       
        function Commit(obj, commit_type)
			% commit time dependent state on progressing to the next time
			% increment
            for m=1:length(obj.models)
                obj.models{m}.Commit(obj, commit_type);
            end
            
            if (commit_type == "Timedep")
                obj.StateVec_Old = obj.StateVec;
				for i=1:obj.dofSpace.NSteps
					obj.K{i} = [];
				end
            end
        end
        
        function anyIrr = Irreversibles(obj)
			% Checks if any irreversable processes can occur at the end of
			% the time step 
            anyIrr = false;
            for m=1:length(obj.models)
                anyIrr = anyIrr + obj.models{m}.Irreversibles(obj);
            end
        end
            
        function Constrain(obj, stp)
			%constrain the tangential stiffness matrix and force vector,
			%based on defined constraints

            obj.convals_corr = obj.convals{stp} - obj.StateVec{stp}(obj.condofs{stp});
            basemat = speye(size(obj.K{stp}));
            obj.unconMat = basemat;
            obj.unconMat(:, obj.condofs{stp}) = [];
            obj.conMat = basemat(:, obj.condofs{stp});

			if (isempty(obj.convals_corr))

			else
            	obj.fint{stp} = obj.unconMat'*obj.fint{stp} + obj.unconMat'*obj.K{stp}*obj.conMat*obj.convals_corr;
            	obj.K{stp}    = obj.unconMat'*obj.K{stp}*obj.unconMat;
			end
        end
        
        function Update(obj, dx, stp)
			% adds increment dx to the state vector
			if (isempty(obj.convals_corr))
				obj.StateVec{stp} = obj.StateVec{stp} + obj.unconMat*dx;
			else
				obj.StateVec{stp} = obj.StateVec{stp} + obj.unconMat*dx + obj.conMat*obj.convals_corr;
			end
            
        end
        
        function info = Request_Info(obj, var, elems, loc)
			%inter-model communicator able to request information based on
			%variable name var
			
            info = false;
            for m=1:length(obj.models)
                [hasInfo, provided] = obj.models{m}.Provide_Info(obj, var, elems, loc);
                if (hasInfo)
                    info = provided;
                end
            end   
        end

    end
end

