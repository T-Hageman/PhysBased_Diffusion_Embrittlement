function [n, t] =  getNormals(obj, group, elem)
%GETNORMALS gets normal vectors in integration points
%   Detailed explanation goes here

    myNode = obj.Elementgroups{group}.Elems(elem,:);
    X = obj.Nodes(myNode,1);
    Y = obj.Nodes(myNode,2);
    
    [n, t] = obj.Elementgroups{group}.ShapeFunc.getNormals(X,Y);

end

