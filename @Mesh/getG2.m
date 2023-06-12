function [G2] = getG2(obj, group, elem)
    % returns shape functions for the element in "group" with number "elem"
    % in the format: N(ip, shape), G(ip, shape, dx/dy), w(ip)
    
    
    myNode = obj.Elementgroups{group}.Elems(elem,:);
    X = obj.Nodes(myNode,1);
    Y = obj.Nodes(myNode,2);
    
    G2 = obj.Elementgroups{group}.ShapeFunc.getG2(X,Y);

end

