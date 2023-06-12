function Propagate_Disc(obj, physics, tipnode, direction)

    neighbours = obj.getConnected("Internal", tipnode);

    if (abs(direction(1))<1e-6 && abs(direction(2)--1)<1e-6) % propagate downwards
        elem_nodes = obj.Elementgroups{1}.Elems(neighbours(1),:);
        SideNodes = elem_nodes([9 6 3]);

        newnodes = substitude(obj, SideNodes(1:2), "Vertical");

        new_nodes_int = [SideNodes(1:3),newnodes,SideNodes(3)];
        obj.Elementgroups{6}.Elems(end+1,:) = new_nodes_int';

        dofTypeIndices = physics.dofSpace.getDofType({"dx","dy","pd"})  ;  

        physics.dofSpace.addDofs(dofTypeIndices(1:2), newnodes);
        for i=1:2
            physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodes(1:2)));
            physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodes(1:2)));  
        end

        physics.dofSpace.addDofs(dofTypeIndices(3), SideNodes(2:3) );
        physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodes(2:3) )) = 0.0;
        physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodes(2:3) )) = 0.0;
    elseif (abs(direction(1)--1)<1e-6 && abs(direction(2))<1e-6) % propagate left
        elem_nodes = obj.Elementgroups{1}.Elems(neighbours(1),:);
        if (size(neighbours,1) == 2)
            SideNodes = elem_nodes([7 8 9]);
        else
            SideNodes = elem_nodes([7 8 9]);
        end
        
        
        newnodes = substitude(obj, SideNodes(2:3), "Horizontal");
        
        new_nodes_int = [SideNodes(1:3),SideNodes(1),newnodes];
        obj.Elementgroups{6}.Elems(end+1,:) = new_nodes_int';
        
        dofTypeIndices = physics.dofSpace.getDofType({"dx","dy","pd"})  ;  
        physics.dofSpace.addDofs(dofTypeIndices(1:2), newnodes);
        for i=1:2
            physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodes(2:3)));
            physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodes(2:3)));  
        end
        
        physics.dofSpace.addDofs(dofTypeIndices(3), SideNodes(1:2) );
        physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodes(1:2) )) = 0.0;
        physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodes(1:2) )) = 0.0;
        
    elseif (abs(direction(1)-1)<1e-6 && abs(direction(2))<1e-6) % propagate right
        elem_nodes = obj.Elementgroups{1}.Elems(neighbours(3),:);
        SideNodes = elem_nodes([7 8 9]);
        
        newnodes = substitude(obj, SideNodes(1:2), "Horizontal");
        new_nodes_int = [SideNodes(1:3),newnodes,SideNodes(3)];
        obj.Elementgroups{6}.Elems(end+1,:) = new_nodes_int';
        
        dofTypeIndices = physics.dofSpace.getDofType({"dx","dy","pd"})  ;  
        physics.dofSpace.addDofs(dofTypeIndices(1:2), newnodes);
        for i=1:2
            physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodes(1:2)));
            physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodes(1:2)));  
        end
        
        physics.dofSpace.addDofs(dofTypeIndices(3), SideNodes(2:3) );
        physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodes(2:3) )) = 0.0;
        physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodes(2:3) )) = 0.0;
     elseif (abs(direction(1)-42)<1e-6 && abs(direction(2)--42)<1e-6) % T-junction
        %split tip 
        newnode = substitude(obj, tipnode, "Tip");
        dofTypeIndices = physics.dofSpace.getDofType({"dx","dy","pd"})  ;  
        physics.dofSpace.addDofs(dofTypeIndices(1:2), newnode);
        for i=1:2
            physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnode)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), tipnode));
            physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnode)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), tipnode));  
        end        
        
        
        %new elements
        elem_nodes  = obj.Elementgroups{1}.Elems(neighbours(1),:);
        elem_nodes2 = obj.Elementgroups{1}.Elems(neighbours(3),:);
        
        SideNodesLeft  = elem_nodes([7 8 9]);
        SideNodesRight = elem_nodes2([7 8 9]);
        
        newnodesLeft = substitude(obj, SideNodesLeft(2:3), "Horizontal");
        new_nodes_intLeft = [SideNodesLeft(1:3),SideNodesLeft(1), newnodesLeft];
        obj.Elementgroups{6}.Elems(end+1,:) = new_nodes_intLeft';
        
        newnodesRight = substitude(obj, SideNodesRight(2), "Horizontal");
        new_nodes_intRight = [SideNodesRight(1:3), newnode, newnodesRight, SideNodesRight(3)];
        obj.Elementgroups{6}.Elems(end+1,:) = new_nodes_intRight';
        
        physics.dofSpace.addDofs(dofTypeIndices(1:2), newnodesLeft);
        for i=1:2
            physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodesLeft)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodesLeft(2:3)));
            physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodesLeft)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodesLeft(2:3)));  
        end  
        
        physics.dofSpace.addDofs(dofTypeIndices(1:2), newnodesRight);
        for i=1:2
            physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodesRight)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodesRight(2:2)));
            physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodesRight)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodesRight(2:2)));  
        end
        
        
        physics.dofSpace.addDofs(dofTypeIndices(3), SideNodesLeft(1:3) );
        physics.dofSpace.addDofs(dofTypeIndices(3), SideNodesRight(2:3) );
        
        physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodesLeft(1:3) )) = 0.0;
        physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodesLeft(1:3) )) = 0.0;
        
        physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodesRight(2:3) )) = 0.0;
        physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodesRight(2:3) )) = 0.0;
        
        % last corrections
        obj.Tjunction = [new_nodes_intLeft(3) new_nodes_intLeft(6)];
        obj.Elementgroups{6}.Elems(end-2,3) = new_nodes_intLeft(6);
        
        physics.dofSpace.addDofs(dofTypeIndices(3), new_nodes_intLeft(6) );
        physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), new_nodes_intLeft(6) )) = 0.0;
        physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), new_nodes_intLeft(6) )) = 0.0;
    else
        
       error("still to implement direction for propagation"); 
        
    end


    obj.check();
end

function newnodes = substitude(obj, oldnodes, direction)
    ToSplit = oldnodes;
    for j=1:length(ToSplit)
        nnodes = size(obj.Nodes, 1) +1;
        obj.Nodes(nnodes,:) = obj.Nodes(ToSplit(j),:);
        for g=1:6
            for el=1:size(obj.Elementgroups{g}.Elems, 1)
                nodeloc = find(obj.Elementgroups{g}.Elems(el,:)==ToSplit(j));
                if (length(nodeloc) >0)
                    for k=1:length(nodeloc)
                        if (direction == "Vertical")
                            if (obj.Elementgroups{g}.type == "Q9" && (nodeloc(k) == 1 || nodeloc(k) == 4 || nodeloc(k) == 7))
                                obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
                            end
                            if (obj.Elementgroups{g}.type == "L3" && nodeloc(k) == 3 )
                                obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
                            end
                            if (obj.Elementgroups{g}.type == "LI6" && nodeloc(k) == 6)
                                obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
                            end
                        elseif (direction == "Horizontal")
                            if (obj.Elementgroups{g}.type == "Q9" && (nodeloc(k) == 1 || nodeloc(k) == 2 || nodeloc(k) == 3))
                                obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
                            end
                            if (obj.Elementgroups{g}.type == "L3" && (nodeloc(k) == 3 || nodeloc(k) == 1))
                                obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
                            end
                            if (obj.Elementgroups{g}.type == "LI6" && (nodeloc(k) == 6 || nodeloc(k) == 4))
                                obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
                            end
                        elseif (direction == "Tip")
                            if (obj.Elementgroups{g}.type == "Q9" && nodeloc(k) == 1)
                                obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
                            end
                            if (obj.Elementgroups{g}.type == "LI6" && (nodeloc(k) == 6 ))
                                obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
                            end
                        else
                           error("Fracture propagation error") 
                        end
                    end
                end
            end
        end
        newnodes(j) = nnodes;
    end

end