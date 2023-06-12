addpath(genpath('./Models'))
addpath(genpath('./Shapes'))

Folders = {"./Results/Subgrid_l_0.000125/";
		 		"./Results/Subgrid_l_0.00025/";
		 		"./Results/Subgrid_l_0.000375/";
		 		"./Results/Subgrid_l_0.0005/";
 		 		"./Results/WuLorenzis_l_0.000125/";
		 		"./Results//WuLorenzis_l_0.00025/";
 		 		"./Results//WuLorenzis_l_0.000375/";
 		 		"./Results//WuLorenzis_l_0.0005/";};

nmax = 440;

sname = {"p125"; "p250"; "p375"; "p500"; "d125"; "d250"; "d375"; "d500"};
%try
	for f=1:length(Folders)
		vid = VideoWriter("Animations/"+sname{f}+".mp4");
		%vid2= VideoWriter("Animations/"+sname{f}+"_frac.avi");
		vid.FrameRate = 10;
		%vid2.FrameRate= 2;
		open(vid);
		%open(vid2);
	
		for i=2:2:nmax
			i
			load(Folders{f}+string(i)+".mat")	
    		figure(42)
			clf 
    		tiledlayout('flow');
			nexttile
        		physics.PlotNodal("phi",0, "Internal");
        		title("\phi")
				colorbar
				axis image
				axis off
			nexttile
        		physics.PlotNodal("CL",0, "Internal");
        		title("C_L")
				colorbar
				axis image
				axis off

			plotElectrolyte(physics)

			set(gcf, 'WindowState','maximized')
			drawnow()
			frame = getframe(gcf);
			writeVideo(vid,frame);
	
	% 		figure(43)
	% 		clf
	% 		physics.models{6}.plotFields(physics);
	% 		set(gcf, 'WindowState','maximized')
	% 		frame = getframe(gcf);
	% 		writeVideo(vid2,frame);
		end
	
	
	end
%catch ME

%end
	close(vid);
% 	close(vid2);


function plotElectrolyte(physics)
	edges = [0 10 10 0 0;
		     0 0 10 10 0]*1e-3;

	nexttile
		physics.models{6}.plotFieldspH(physics);
		title("pH")
		colorbar
		axis image
		axis off
		hold on
		plot(edges(1,:),edges(2,:),'k')

	nexttile
		physics.models{6}.plotFieldsFe(physics);
		title("C_{Fe^{2+}}")
		colorbar
		axis image
		axis off
		hold on
		plot(edges(1,:),edges(2,:),'k')


end

