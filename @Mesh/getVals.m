function [N, G, w] = getVals(obj, group, elem)
    % returns shape functions for the element in "group" with number "elem"
    % in the format: N(ip, shape), G(ip, shape, dx/dy), w(ip)
    
    
    myNode = obj.Elementgroups{group}.Elems(elem,:);
    X = obj.Nodes(myNode,1);
    Y = obj.Nodes(myNode,2);
    
    [N, G, w] = obj.Elementgroups{group}.ShapeFunc.getVals(X,Y);

end

