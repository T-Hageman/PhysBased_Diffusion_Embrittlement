function [G2] = getG2(obj, group, elem)
    %GETG2 returns second gradient of shape functions for the element in "group" with number "elem"
    % in the format: G2(ip, shape, dxx/dyy/dxy)
    
    
    myNode = obj.Elementgroups{group}.Elems(elem,:);
    X = obj.Nodes(myNode,1);
    Y = obj.Nodes(myNode,2);
    
    G2 = obj.Elementgroups{group}.ShapeFunc.getG2(X,Y);

end

