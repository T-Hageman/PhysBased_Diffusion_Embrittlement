classdef DummyPhaseField < BaseModel
    %DummyPhaseField Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        l
		Phi_Step
		dx_Step
    end
    
    methods
        function obj = DummyPhaseField(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "DummyPhaseField";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
			[obj.dofTypeIndices, stp]  = obj.dofSpace.getDofType({"dx","dy","phi"});
			obj.dx_Step = stp(1);
			obj.Phi_Step= stp(3);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            %% get parameters
            obj.l = inputs.l;
        end
        
        function getKf(obj, physics, stp)

			if (stp == obj.Phi_Step)
				Hdom = 10e-3;
	
            	fprintf("        DummyPhaseField get Matrix:")
            	t = tic;
            	
            	dt = physics.dt;
	
            	allNodes = obj.mesh.GetAllNodesForGroup(obj.myGroupIndex);
				newcons = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), allNodes);
            	nodeCons = 0*allNodes;
            	for i=1:length(allNodes)
					xy = [obj.mesh.Nodes(allNodes(i),1), obj.mesh.Nodes(allNodes(i),2)];
					
					if (xy(1)<5e-3)
						dst = Hdom/2-xy(2);
					else
						dst=sqrt((xy(1)-Hdom/2)^2+(xy(2)-5e-3)^2);
					end
					
					if (abs(dst)>5.001*obj.l)
						pf = 0;
					else
						pf = exp(-abs(dst)/(2*obj.l));
					end
					nodeCons(i) = pf;
	
	
            	end
            	
				physics.condofs{obj.Phi_Step} = [physics.condofs{obj.Phi_Step}; newcons];
            	physics.convals{obj.Phi_Step} = [physics.convals{obj.Phi_Step}; nodeCons];
				physics.StateVec{obj.Phi_Step}(newcons) = nodeCons;
				physics.StateVec_Old{obj.Phi_Step}(newcons) = nodeCons;
	
% 				figure(751)
% 				clf
% 				physics.PlotNodal("phi",1000, "Internal");
            	
            	tElapsed = toc(t);
            	fprintf("            (Assemble time:"+string(tElapsed)+")\n");
			end
        end

        function [hasInfo, provided] = Provide_Info(obj, physics, var, elems, loc)
           hasInfo = false;
           provided = [];
           
        end
        
    end
end

