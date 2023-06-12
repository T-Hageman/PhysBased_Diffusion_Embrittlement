classdef Constrainer < BaseModel
    %CONSTRAINER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        myName
        mesh
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        conVal
		conStep
    end
    
    methods
        function obj = Constrainer(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "Constrainer";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Ngroup;
            obj.myGroupIndex = obj.mesh.getNodeGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            [obj.dofTypeIndices, obj.conStep] = obj.dofSpace.getDofType(inputs.dofs);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.Nodegroups{obj.myGroupIndex}.Nodes);
            
            obj.conVal = inputs.conVal;
        end
        
        function getKf(obj, physics, stp)
            allNodes = obj.mesh.Nodegroups{obj.myGroupIndex}.Nodes;
            
            for i=1:length(obj.dofTypeIndices)
				if (obj.conStep == stp)
                	newcons = obj.dofSpace.getDofIndices(obj.dofTypeIndices(i), allNodes);
                	
                	physics.condofs{stp} = [physics.condofs{stp}; newcons];
                	physics.convals{stp} = [physics.convals{stp}; newcons*0+obj.conVal(i)];
				end
            end
        end
    end
end

