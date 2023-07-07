classdef DofSpace < handle
    %DOFSPACE Management of degrees of freedom
    %   Class which keeps track of added degrees of freedom, giving index
	%   of dof for a dof name, and location for a combination of index and
	%   node.
    
    properties
        mesh			%pointer to mesh object
        DofTypes		%String vector containing names of each dof type
        DofSteps		%int vector containing the steps in which each dof type is solved for
        DofNumbering	%node, dofType -> dofnumber
        NDofs			%total amout of degrees of freedom
		NSteps			%number of solution steps
    end
    
    methods
        function obj = DofSpace(mesh, dofs_in)
            %DOFSPACE Construct an instance of this class
            
            obj.NSteps = max(dofs_in.Step);

            obj.mesh = mesh;
            obj.DofTypes = dofs_in.dofs;
            obj.DofSteps = dofs_in.Step;
            obj.NDofs = zeros(obj.NSteps,1);
            obj.DofNumbering = sparse(length(obj.mesh.Nodes), length(obj.DofTypes));
        end

        function addDofs(obj, dofIndices, nodeIndex)
			% Adds degrees of freedom for type "dofIndices" to the nodes
			% "nodeIndex"

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
			% returns the dof type index for pre-existing degrees of
			% freedom

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
           % gets the indices for a combination of degree of freedom
		   % "doftype" and nodes "NodeIndices"  
		   
			DofIndices = obj.DofNumbering(NodeIndices, dofType);
			if (isempty(DofIndices==0) && length(DofIndices)>=1)
				disp("error here (dofspace)") 
			end
		end

    end
end

