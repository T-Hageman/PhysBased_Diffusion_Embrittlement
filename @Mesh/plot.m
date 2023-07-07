function f1 = plot(obj, plotnodes, plotelems, plotnames)
	%PLOT Plot the used mesh to a new window, and show node and element
	%numbering (extremely slow, use for debugging only)

    f1 = figure();
    if (plotnodes)
        plot(obj.Nodes(:,1),obj.Nodes(:,2), 'k*');
    end 
    hold on
    if (plotelems)
        for g=1:length(obj.Elementgroups)
			clear X Y textloc textlabel
            if (obj.Elementgroups{g}.type == "Q9")
                for el=1:size(obj.Elementgroups{g}.Elems, 1)
                    elnodes = obj.Elementgroups{g}.Elems(el,:);
                    order = [1 2 3 6 9 8 7 4];
                    X(el,:) = obj.Nodes(elnodes(order),1);
                    Y(el,:) = obj.Nodes(elnodes(order),2);

                    textloc(el,:) = [obj.Nodes(elnodes(5),1) obj.Nodes(elnodes(5),2)];
                    textlabel(el) = "E"+string(el);
                end
                patch('XData',X','YData',Y','EdgeColor','black','FaceColor','none');
                hold on
                if (plotnames)
                   text(textloc(:,1), textloc(:,2), textlabel);
                end
			end
			if (obj.Elementgroups{g}.type == "T6")
                for el=1:size(obj.Elementgroups{g}.Elems, 1)
                    elnodes = obj.Elementgroups{g}.Elems(el,:);
                    order = [1 4 2 5 3 6];
                    X(el,:) = obj.Nodes(elnodes(order),1);
                    Y(el,:) = obj.Nodes(elnodes(order),2);

                    textloc(el,:) = [(obj.Nodes(elnodes(1),1)+obj.Nodes(elnodes(2),1)+obj.Nodes(elnodes(3),1))/3 ...
						             (obj.Nodes(elnodes(1),2)+obj.Nodes(elnodes(2),2)+obj.Nodes(elnodes(3),2))/3];
                    textlabel(el) = "E"+string(el);
                end
                patch('XData',X','YData',Y','EdgeColor','black','FaceColor','none');
                hold on
                if (plotnames)
                   text(textloc(:,1), textloc(:,2), textlabel);
                end
            end
            clear X Y order
            if (obj.Elementgroups{g}.type == "L3")
                for el=1:size(obj.Elementgroups{g}.Elems, 1)
                    elnodes = obj.Elementgroups{g}.Elems(el,:);
                    order = [1 2 3];
                    X(el,:) = obj.Nodes(elnodes(order),1);
                    Y(el,:) = obj.Nodes(elnodes(order),2);
                end
                plot(X,Y,'r')
                hold on
            end
            if (obj.Elementgroups{g}.type == "LI6")
                for el=1:size(obj.Elementgroups{g}.Elems, 1)
                    elnodes = obj.Elementgroups{g}.Elems(el,:);
                    order = [1 2 3 6 5 4];
                    X(el,:) = obj.Nodes(elnodes(order),1);
                    Y(el,:) = obj.Nodes(elnodes(order),2);
                end
                plot(X,Y,'b')
                hold on
            end
        end
    end
end

