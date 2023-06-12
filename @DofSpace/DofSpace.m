classdef DofSpace < handle
    %DOFSPACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        DofTypes
        DofSteps
        DofNumbering	%node, dofType
        NDofs			%step
		NSteps
    end
    
    methods
        function obj = DofSpace(mesh, dofs_in)
            %DOFSPACE Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.NSteps = max(dofs_in.Step);

            obj.mesh = mesh;
            obj.DofTypes = dofs_in.dofs;
            obj.DofSteps = dofs_in.Step;
            obj.NDofs = zeros(obj.NSteps,1);
            obj.DofNumbering = sparse(length(obj.mesh.Nodes), length(obj.DofTypes));
        end

        function addDofs(obj, dofIndices, nodeIndex)

            for i=1:length(dofIndices)
                dofStep = obj.DofSteps(dofIndices(i));
				
                for k=1:length(nodeIndex)
                    if (nodeIndex(k)<=size(obj.DofNumbering, 1))
                        curNum = obj.DofNumbering(nodeIndex(k), dofIndices(i));
                    else
                        obj.DofNumbering(nodeIndex(k), :) = 0;
                        curNum = 0;
                    end
                    if (curNum == 0) %dof does not yet exist if 0, add dof
                        obj.NDofs(dofStep) = obj.NDofs(dofStep)+1;
                        obj.DofNumbering(nodeIndex(k), dofIndices(i)) = obj.NDofs(dofStep);
                    end
                end
            end
        end
        
		function [DofTypeIndex, DofStepIndex] = getDofType(obj, dofnames)
            DofTypeIndex = zeros(length(dofnames),1);
			DofStepIndex = zeros(length(dofnames),1);

            for i=1:length(dofnames)
               for j=1:length(obj.DofTypes)
                  if (dofnames{i}==obj.DofTypes{j})
                      DofTypeIndex(i) = j;
					  DofStepIndex(i) = obj.DofSteps(j);
                  end
               end
            end
        end
        
		function DofIndices = getDofIndices(obj, dofType, NodeIndices)
           
			DofIndices = obj.DofNumbering(NodeIndices, dofType);
			if (isempty(DofIndices==0) && length(DofIndices)>=1)
				disp("error here (dofspace)") 
			end
		end

    end
end

