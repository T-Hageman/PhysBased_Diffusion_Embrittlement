function xy = getIPCoords(obj, group, elem)
    %GETIPCOORDS returns integration point coordinates for the element in "group" with number "elem"
        
    myNode = obj.Elementgroups{group}.Elems(elem,:);
    X = obj.Nodes(myNode,1);
    Y = obj.Nodes(myNode,2);
    
    xy = obj.Elementgroups{group}.ShapeFunc.getIPGlobal(X,Y);
end

