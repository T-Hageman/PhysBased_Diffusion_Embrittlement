function PlotIP(obj, varName, plotloc)
	%PLOTIP Plots integration-point level data
	
    for g=1:length(obj.mesh.Elementgroups)
        if (obj.mesh.Elementgroups{g}.name == plotloc)
            for e=1:size(obj.mesh.Elementgroups{g}.Elems, 1)
                IPvar = obj.Request_Info(varName, e, "Interior");
                xy = obj.mesh.getIPCoords(g, e);
                
				meanx = mean(xy(1,:));
				meany = mean(xy(2,:));
				for i=1:size(xy,2)
					if xy(1,i)<meanx
						xy(1,i) = xy(1,i) + 1e-9;
					else
						xy(1,i) = xy(1,i) - 1e-9;
					end
					if xy(2,i)<meany
						xy(2,i) = xy(2,i) + 1e-9;
					else
						xy(2,i) = xy(2,i) - 1e-9;
					end
				end

                x(e, :) = xy(1,:);
                y(e, :) = xy(2,:);
                z(e, :) = IPvar;
            end
            
            xlims = [min(min(x)), max(max(x))];
            ylims = [min(min(y)), max(max(y))];
            
            [xi,yi] = meshgrid(linspace(xlims(1), xlims(2), 1000), linspace(ylims(1), ylims(2), 1000));
            zi = griddata(x,y,z,xi,yi,'linear');
            %s = pcolor(xi,yi,zi);
			s = surf(xi, yi, zi);
            s.EdgeColor = 'none';
            hold on
            colorbar
        end
    end
end


