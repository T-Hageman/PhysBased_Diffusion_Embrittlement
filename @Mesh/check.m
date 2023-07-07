function check(obj)
	%CHECK Calculate the area of the mesh and print it to confirm shape
	%functions for all elements can be evaluated, and (manually) confirm
	%that the total mesh area is correct

    disp("Checking mesh");

    newArea = 0.0*obj.Area;
    for g=1:length(obj.Elementgroups)
        for ielem = 1:size(obj.Elementgroups{g}.Elems, 1)
            [N, G, w] = obj.getVals(g, ielem);
            newArea(g) = newArea(g) + sum(w);
        end
        
        disp("    "+obj.Elementgroups{g}.name+ ": Old Area: "+string(obj.Area(g))+"    new Area: "+string(newArea(g)))
    end


    obj.Area = newArea;
end

