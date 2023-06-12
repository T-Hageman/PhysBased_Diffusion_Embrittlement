function [Nodes, Elementgroups, Nodegroups, Area, rect] = Square_Generator_WithDisc(obj, props)
    rect = true;

    Lx = props.Lx; 
    Ly = props.Ly;  
    Nx = props.Nx;  
    Ny = props.Ny;  
	HFrac = props.Uextern;

	%bottom half
    nnodes = 0;
    nnodes_per_row = Nx*2+1;
    for j=0:Ny
        y = Ly*j/(Ny*2);
        for i=0:Nx*2
            x = Lx*i/(Nx*2);
            
            nnodes = nnodes+1;
            Nodes(nnodes, 1:2) = [x,y];
        end
    end

    Elementgroups{1}.name = "Metal";
    Elementgroups{1}.type = "Q9";
    Area(1) = Lx*Ly;
    
    cntr = 0;
    for j=1:Ny/2
        for i=1:Nx
            row1 = (2*(j-1)*nnodes_per_row)+(2*(i-1))+[1 2 3];
            elem_nodes = [row1 row1+nnodes_per_row row1+2*nnodes_per_row];
            
            cntr = cntr+1;
            Elementgroups{1}.Elems(cntr,:) = elem_nodes;
        end
    end
    
    Elementgroups{2}.name = "M_Bottom";
    Elementgroups{2}.type = "L3";
    Area(2) = Lx;
    
    Elementgroups{3}.name = "Interface";
    Elementgroups{3}.type = "L3";
    Area(3) = Lx;
    
    cntr = 0;
	for i=1:Nx
        top_nodes = ((Ny)*nnodes_per_row)+(2*(i-1))+[1 2 3];
        bot_nodes = (2*(i-1))+[1 2 3];
            
        cntr = cntr+1;
        Elementgroups{2}.Elems(cntr,:) = bot_nodes;
        Elementgroups{3}.Elems(cntr,:) = top_nodes(length(bot_nodes):-1:1);
    end
    
    Elementgroups{4}.name = "M_Left";
    Elementgroups{4}.type = "L3";
    Area(4) = Ly;
    
    Elementgroups{5}.name = "M_Right";
    Elementgroups{5}.type = "L3";
    Area(5) = Ly;
    
    cntr = 0;
	for i=1:Ny/2
        left_nodes =  (2*(i-1))*nnodes_per_row  +   [1 nnodes_per_row+1 2*nnodes_per_row+1];
        right_nodes = (2*(i-1))*nnodes_per_row  +   [nnodes_per_row 2*nnodes_per_row 3*nnodes_per_row];
            
        cntr = cntr+1;
        Elementgroups{4}.Elems(cntr,:) = left_nodes(length(left_nodes):-1:1);
        Elementgroups{5}.Elems(cntr,:) = right_nodes;
	end

	%Electrolyte Layer
	NyFrac = 10;
    for j=1:2*(NyFrac)
        y = Ly/2+HFrac*j/(2*NyFrac);
        for i=0:Nx*2
            x = Lx*i/(Nx*2);
            
            nnodes = nnodes+1;
            Nodes(nnodes, 1:2) = [x,y];
        end
    end

    Elementgroups{6}.name = "Electrolyte";
    Elementgroups{6}.type = "Q9";
    Area(6) = Lx*HFrac;

    cntr = 0;
	Elems_Left = [];
    for j=Ny/2:Ny/2+NyFrac-1
        for i=1:Nx
            row1 = ((2*(j))*nnodes_per_row)+(2*(i-1))+[1 2 3];
            elem_nodes = [row1 row1+nnodes_per_row row1+2*nnodes_per_row];
            
            cntr = cntr+1;
            Elementgroups{6}.Elems(cntr,:) = elem_nodes;
			if (i==1) 
				Elems_Left = [Elems_Left; elem_nodes(1), elem_nodes(4), elem_nodes(7)]; 
			end
        end
    end

    Elementgroups{7}.name = "Interface2";
    Elementgroups{7}.type = "L3";
    Area(7) = Lx;
    
    cntr = 0;
	for i=1:Nx
        top_nodes = ((nnodes/nnodes_per_row-1)*nnodes_per_row)+(2*(i-1))+[1 2 3];
            
        cntr = cntr+1;
        Elementgroups{7}.Elems(cntr,:) = top_nodes(length(bot_nodes):-1:1);
	end

    Elementgroups{8}.name = "E_Left";
    Elementgroups{8}.type = "L3";
	Area(8) = HFrac;
	Elementgroups{8}.Elems = Elems_Left;

	
	%top half
    for j=1:Ny
        y = Ly/2+HFrac+Ly*j/(Ny*2);
        for i=0:Nx*2
            x = Lx*i/(Nx*2);
            
            nnodes = nnodes+1;
            Nodes(nnodes, 1:2) = [x,y];
        end
	end
    
    cntr = length(Elementgroups{1}.Elems);
    for j=Ny/2+NyFrac+1:Ny+NyFrac
        for i=1:Nx
            row1 = ((2*(j-1))*nnodes_per_row)+(2*(i-1))+[1 2 3];
            elem_nodes = [row1 row1+nnodes_per_row row1+2*nnodes_per_row];
            
            cntr = cntr+1;
            Elementgroups{1}.Elems(cntr,:) = elem_nodes;
        end
	end

	Elementgroups{9}.name = "M_Top";
    Elementgroups{9}.type = "L3";
    Area(9) = Lx;
    
    cntr = 0;
	for i=1:Nx
        top_nodes = (2*(Ny+NyFrac)*nnodes_per_row)+(2*(i-1))+[1 2 3];
            
        cntr = cntr+1;
        Elementgroups{9}.Elems(cntr,:) = top_nodes(length(bot_nodes):-1:1);
	end

	%nodegroups
	for g=1:length(Elementgroups)
        Nodegroups{g}.name = Elementgroups{g}.name;
        Nodegroups{g}.Nodes = unique(reshape(Elementgroups{g}.Elems,[],1));
	end

end

