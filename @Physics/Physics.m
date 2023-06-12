classdef Physics < handle
    %PHYSICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        models
        dofSpace
        
        K
        fint
        StateVec
        StateVec_Old
        time
        dt
        
        condofs
        convals
        convals_corr
        conMat
        unconMat
		nonz
    end
    
    methods
        PlotNodal(obj, dofName, dispscale, plotloc)
        PlotIP(obj, varName, plotloc)
        
        function obj = Physics(mesh, inputs, dofs_in)
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
            for m=1:length(obj.models)
                obj.models{m}.OncePerStep(obj, stp);
			end
		end


        function Assemble(obj, stp)
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
            anyIrr = false;
            for m=1:length(obj.models)
                anyIrr = anyIrr + obj.models{m}.Irreversibles(obj);
            end
        end
            
        function Constrain(obj, stp)
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
			if (isempty(obj.convals_corr))
				obj.StateVec{stp} = obj.StateVec{stp} + obj.unconMat*dx;
			else
				obj.StateVec{stp} = obj.StateVec{stp} + obj.unconMat*dx + obj.conMat*obj.convals_corr;
			end
            
        end
        
        function info = Request_Info(obj, var, elems, loc)
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

