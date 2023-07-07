function [Nodes, Elementgroups, Nodegroups, Area, rect] = Square_Generator(obj, props)
    rect = true;

    Lx = props.Lx; 
    Ly = props.Ly;  
    Nx = props.Nx;  
    Ny = props.Ny;  

    nnodes = 0;
    nnodes_per_row = Nx*2+1;
    for j=0:Ny*2
        y = Ly*j/(Ny*2);
        for i=0:Nx*2
            x = Lx*i/(Nx*2);
            
            nnodes = nnodes+1;
            Nodes(nnodes, 1:2) = [x,y];
        end
    end


    Elementgroups{1}.name = "Internal";
    Elementgroups{1}.type = "Q9";
    Area(1) = Lx*Ly;
    
    cntr = 0;
    for j=1:Ny
        for i=1:Nx
            row1 = (2*(j-1)*nnodes_per_row)+(2*(i-1))+[1 2 3];
            elem_nodes = [row1 row1+nnodes_per_row row1+2*nnodes_per_row];
            
            cntr = cntr+1;
            Elementgroups{1}.Elems(cntr,:) = elem_nodes;
        end
    end
    
    Elementgroups{2}.name = "Bottom";
    Elementgroups{2}.type = "L3";
    Area(2) = Lx;
    
    Elementgroups{3}.name = "Top";
    Elementgroups{3}.type = "L3";
    Area(3) = Lx;
    
    cntr = 0;
	for i=1:Nx
        top_nodes = ((2*Ny)*nnodes_per_row)+(2*(i-1))+[1 2 3];
        bot_nodes = (2*(i-1))+[1 2 3];
            
        cntr = cntr+1;
        Elementgroups{2}.Elems(cntr,:) = bot_nodes;
        Elementgroups{3}.Elems(cntr,:) = top_nodes(length(bot_nodes):-1:1);
    end
    
    Elementgroups{4}.name = "Left";
    Elementgroups{4}.type = "L3";
    Area(4) = Ly;
    
    Elementgroups{5}.name = "Right";
    Elementgroups{5}.type = "L3";
    Area(5) = Ly;
    
    cntr = 0;
	for i=1:Ny
        left_nodes =  (2*(i-1))*nnodes_per_row  +   [1 nnodes_per_row+1 2*nnodes_per_row+1];
        right_nodes = (2*(i-1))*nnodes_per_row  +   [nnodes_per_row 2*nnodes_per_row 3*nnodes_per_row];
            
        cntr = cntr+1;
        Elementgroups{4}.Elems(cntr,:) = left_nodes(length(left_nodes):-1:1);
        Elementgroups{5}.Elems(cntr,:) = right_nodes;
	end
    
    for g=1:length(Elementgroups)
        Nodegroups{g}.name = Elementgroups{g}.name;
        Nodegroups{g}.Nodes = unique(reshape(Elementgroups{g}.Elems,[],1));
	end

    Nodegroups{g+1}.name = "LeftBottom";
    Nodegroups{g+1}.Nodes = 1;

	Nodegroups{g+2}.name = "LeftTop";
    Nodegroups{g+2}.Nodes = min(min(Elementgroups{3}.Elems));
end

