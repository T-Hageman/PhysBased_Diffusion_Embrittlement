classdef Mesh < handle
    %MESH Class that handles the geometric description and discretisation
	% of the domain
    
    properties
        Area			%Vector containing the area/length of each element group
        Elementgroups	%Structure containing all element groups (in turn containing the base element shapes and elements)
        Nodegroups		%Nodegroups
        Nodes			%Vector containing the coordinates for all nodes
        ipcount1D		% number of one-dimensional integration point
        zeroWeight		%flag indicating whether zero-weight integration points have been added
    end
    
    methods
		% functions defined in other files
        f1 = plot(obj, plotnodes, plotelems, plotnames);
        xy = getIPCoords(obj, group, elems);
        check(obj);
        [N, G, w] = getVals(obj, group, elem);
		G2 = getG2(obj, group, elem);
        [n, t] = getNormals(obj, group, elem);
        
        function obj = Mesh(inProps)
             %MESH Construct an instance of this class
            if (inProps.type == "Square")
                [obj.Nodes, obj.Elementgroups, obj.Nodegroups, obj.Area, rectangular] = obj.Square_Generator(inProps);
			end
			if (inProps.type == "Fracture")
				[obj.Nodes, obj.Elementgroups, obj.Nodegroups, obj.Area, rectangular] = obj.Fracture_Generator(inProps);
			end
			if (inProps.type == "Square_WithDisc")
				[obj.Nodes, obj.Elementgroups, obj.Nodegroups, obj.Area, rectangular] = obj.Square_Generator_WithDisc(inProps);
			end
            
            obj.zeroWeight = inProps.zeroWeight;
            obj.ipcount1D = inProps.ipcount1D;
            for g=1:length(obj.Elementgroups)
               f = str2func(obj.Elementgroups{g}.type);
               obj.Elementgroups{g}.ShapeFunc = f(obj.ipcount1D, rectangular, obj.zeroWeight);
            end
            if (obj.zeroWeight)
                obj.ipcount1D = obj.ipcount1D+2;
            end
            
        end
        
        function groupIndex = getGroupIndex(obj, groupname)
			% returns the element group index for a string "groupname"
            groupIndex = -1;
            for i=1:length(obj.Elementgroups)
               if (groupname == obj.Elementgroups{i}.name)
                   groupIndex = i;
               end
            end
            if (groupIndex == -1)
                msg = "Elementgroup "+groupname+" is not defined.";
                error(msg);
            end
		end

		function groupIndex = getNodeGroupIndex(obj, groupname)
			% returns the node group index for a string "groupname"
            groupIndex = -1;
            for i=1:length(obj.Nodegroups)
               if (groupname == obj.Nodegroups{i}.name)
                   groupIndex = i;
               end
            end
            if (groupIndex == -1)
                msg = "Nodegroup "+groupname+" is not defined.";
                error(msg);
            end
        end
        
        function myNode = getNodes(obj, group, elem)
			%Gets the nodes corresponding to the element number elem within
			%group group
            myNode = obj.Elementgroups{group}.Elems(elem,:);
        end
        
        function nds = GetAllNodesForGroup(obj, groupindex)
			%Gets all nodes for an element group
            myElems = obj.Elementgroups{groupindex}.Elems;
            nds = unique(reshape(myElems, [], 1));
        end
        
        
    end
end

