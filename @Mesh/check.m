function check(obj)

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

