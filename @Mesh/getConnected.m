function elems = getConnected(obj, groupname, node)

    groupnum = obj.getGroupIndex(groupname);
    [elem_ind,~] = find(obj.Elementgroups{groupnum}.Elems==node);

    all_nodes = obj.Elementgroups{groupnum}.Elems(elem_ind, :);
    for i=1:size(all_nodes, 1)
        xy = obj.Nodes(all_nodes(i,:),:);
        nodemean(i,:) = mean(xy);
    end
    [~, ind] = sortrows(nodemean);
    
    elems = elem_ind(ind);
    
end

