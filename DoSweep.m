function DoSweep(i_in)
close all

surfplots = true;
figure(1234567)
tiledlayout(3,4);
figure(12345678)
t = tiledlayout(3,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';

	addpath(genpath('./Models'))
	addpath(genpath('./Shapes'))
	
	cspec = {"r","g","b","m"};
	lspec = {"-.","--","-",":"};

	i=0;
	for modeli=1:2    % 2*4*4 = 32x
		m = ["Subgrid";"WuLorenzis"];
		model = m(modeli);
		for ui=1:4
		li = 0;
		for l = 5e-4:2.5e-4:1.25e-3
			li = li+1;

				u=10^(-8+ui);
				
				i=i+1;
				simsets{i}.u=u;
				simsets{i}.l=l;
				simsets{i}.model=model;
				simsets{i}.linespec = cspec{li}+lspec(ui);
			end
		end
	end
	
 	if true
 		for i=i_in:i_in
 			[physics, tvec, CL_vec]  = main_Static(simsets{i}.model,simsets{i}.l,simsets{i}.u);
 		end
 	end
	figure(500002)
	tiledlayout(1,2);
	for i=1:length(simsets)
		sname = simsets{i}.model+'_'+string(simsets{i}.l)+'_'+string(simsets{i}.u);
		savefolder = "./Results_NoProp/"+sname;
		filename = savefolder+"/end";
		dispname = "$\ell="+string(simsets{i}.l*1000)+",\; U_{\mathrm{ext}}="+string(simsets{i}.u*1000)+"$";
		try
    		load(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax","LFrac");
			if (simsets{i}.model == "Subgrid")
				figure(500002)
				nexttile(1);
			else
				figure(500002)
				nexttile(2);
			end
			plot(tvec/3600,CL_vec,simsets{i}.linespec,'DisplayName',dispname,'LineWidth',1.5);
			hold on

			if surfplots
			f1 = figure(i);
				physics.PlotNodal("phi",0, "Internal");
				caxis([0 1]);
				cb = colorbar;
				cb.Title.String = {'$\phi$', '[$-$]'};
				cb.Title.Interpreter='latex';
				axis off
				f1.Units = "centimeters";
				f1.Position = [5,1,3,6];
				cb.Position = [0.8 0.1 0.05 0.7];
				ax = gca;
				ax.Position = [0.05 0.05 0.6 0.9];
			f2 = figure(1234567)
			nexttile
        		physics.PlotNodal("CL",0, "Internal");
				caxis([0 5]);
				%cb = colorbar;
				%cb.Title.String = {'$C_\mathrm{L}$', '[$\mathrm{mol}/\mathrm{m}^3$]'};
				%cb.Title.Interpreter='latex';
				axis off
				%f2.Units = "centimeters";
				%f2.Position = [8,1,3,6];
				%cb.Position = [0.8 0.1 0.05 0.7];
				%				ax = gca;
				%ax.Position = [0.05 0.05 0.6 0.9];
			f3 = figure(12345678);
			nexttile
				physics.models{6}.plotFieldspH(physics);
				caxis([5 6]);
				hold on
				%plot([0 5e-3 5e-3 0 0],[0 0 50e-3 50e-3 0],'k-')
				%cb = colorbar;
				%cb.Title.String = {'$\mathrm{pH}$', '[$-$]'};
				%cb.Title.Interpreter='latex';
				axis off
				%f3.Units = "centimeters";
				%f3.Position = [11,1,3,6];
				%cb.Position = [0.8 0.1 0.05 0.7];
				%				ax = gca;
				%ax.Position = [0.05 0.05 0.6 0.9];

				figure(f1)
					print(gcf, "./Figures/PhaseField_"+sname+".png",'-dpng','-r1200')
					print(gcf, "./Figures/PhaseField_"+sname+".jpg",'-djpeg','-r1200')
    				print(gcf, "./Figures/PhaseField_"+sname+".eps",'-depsc','-r1200')
					print(gcf, "./Figures/PhaseField_"+sname+".emf",'-dmeta','-r1200')
				%figure(f2)
					%print(gcf, "./Figures/CL_"+sname+".png",'-dpng','-r1200')
					%print(gcf, "./Figures/CL_"+sname+".jpg",'-djpeg','-r1200')
    				%print(gcf, "./Figures/CL_"+sname+".eps",'-depsc','-r1200')
					%print(gcf, "./Figures/CL_"+sname+".emf",'-dmeta','-r1200')
				%figure(f3)
					%print(gcf, "./Figures/pH_"+sname+".png",'-dpng','-r1200')
					%print(gcf, "./Figures/pH_"+sname+".jpg",'-djpeg','-r1200')
    				%print(gcf, "./Figures/pH_"+sname+".eps",'-depsc','-r1200')
					%print(gcf, "./Figures/pH_"+sname+".emf",'-dmeta','-r1200')

					close(f1)
					%close(f2)
					%close(f3)
			end
		catch ec
	
		end
	end




	mrkrs = {"kx-","kd-","ko-","ks-","k-"}
	mrkrs2 = {"x","d","o","s"}
	f3 = figure(1);
	subplot(4,1,1);
	for j=1:4
		filename = "./Results_NoProp/"+"Discrete_"+string(10^(-8+j))+"/end";
		load(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax","LFrac");
		figure(500002)
		nexttile(1);
		hp = plot(tvec/3600,CL_vec,mrkrs{j},'DisplayName','Discrete $U_{\mathrm{ext}}='+string(10^(-8+j+3))+"$");
		plotsparsemarkers(hp, [], mrkrs2{j}, 20 ,false)
		nexttile(2);
		hp = plot(tvec/3600,CL_vec,mrkrs{j},'DisplayName','Discrete $U_{\mathrm{ext}}='+string(10^(-8+j+3))+"$");
		plotsparsemarkers(hp, [], mrkrs2{j}, 20 ,false)

		if surfplots
		f2 = figure(1234567);
		nexttile
        	physics.PlotNodal("CL",0, "Metal");	
			caxis([0 5]);
			ylim([0 50e-3]);
			%cb = colorbar;
			%cb.Title.String = {'$C_\mathrm{L}$', '[$\mathrm{mol}/\mathrm{m}^3$]'};
			%cb.Title.Interpreter='latex';
			axis off
			%f2.Units = "centimeters";
			%f2.Position = [8,1,3,6];
			%cb.Position = [0.8 0.1 0.05 0.7];
			%				ax = gca;
			%ax.Position = [0.05 0.05 0.6 0.9];

		figure(12345678);
		nexttile
			physics.models{7}.plotFieldspH(physics);
			caxis([5 6]);
			hold on
			axis off
% 			f3.Units = "centimeters";
% 			f3.Position = [11,1,3,6];
% 			cb.Position = [0.8 0.1 0.05 0.7];
% 							ax = gca;
% 			ax.Position = [0.05 0.45 0.6 0.1];

			sname = "Discrete"+'_'+string(10^(-8+j));
			%figure(f2)
			%	print(gcf, "./Figures/CL_"+sname+".png",'-dpng','-r1200')
			%	print(gcf, "./Figures/CL_"+sname+".jpg",'-djpeg','-r1200')
			%	print(gcf, "./Figures/CL_"+sname+".eps",'-depsc','-r1200')
 			figure(f3)
% 				print(gcf, "./Figures/pH_"+sname+".png",'-dpng','-r1200')
% 				print(gcf, "./Figures/pH_"+sname+".jpg",'-djpeg','-r1200')
% 				print(gcf, "./Figures/pH_"+sname+".eps",'-depsc','-r1200')
		end
	end
	if surfplots
	figure(12345678)
	cb = colorbar;
	cb.Layout.Tile = 'east'; 
	cb.Title.String = {'$\mathrm{pH}$', '[$-$]'};
	cb.Title.Interpreter='latex';
	f3.Units = "centimeters";
	f3.Position = [11,1,14,6];
	cb.Position = [0.8 0.1 0.025 0.75];
	for i=1:4
		subplot(4,1,i)
		ax = gca;
		ax.Position(3) = 0.8*ax.Position(3);
	end

	print(gcf, "./Figures/pH_"+"Discrete"+".png",'-dpng','-r1200')
 	print(gcf, "./Figures/pH_"+"Discrete"+".jpg",'-djpeg','-r1200')
 	print(gcf, "./Figures/pH_"+"Discrete"+".eps",'-depsc','-r1200')

	figure(f2)
		cbh = colorbar(); 
		cbh.Layout.Tile = 'east'; 
		f2.Units = "centimeters"
		f2.Position = [1,1,14,20];

		print(gcf, "./Figures/CL_"+".png",'-dpng','-r1200')
		print(gcf, "./Figures/CL_"+".jpg",'-djpeg','-r1200')
    	print(gcf, "./Figures/CL_"+".eps",'-depsc','-r1200')
		print(gcf, "./Figures/CL_"+".emf",'-dmeta','-r1200')
	end

	f = figure(500002);
	nexttile(1)
	ax1 = gca;
	xlim([0 200]);
	ylim([0 1.5]);
	xlabel("time [hours]",'Interpreter','latex');
	ylabel("$\overline{C}_\mathrm{L}\;[\mathrm{mol}/\mathrm{m}^3]$",'Interpreter','latex')
	nexttile(2)
	ax2 = gca;
	xlim([0 200]);
	ylim([0 1.5]);
	xlabel("$\mathrm{time}\;[\mathrm{hours}]$",'Interpreter','latex');
	ylabel("$\overline{C}_\mathrm{L}\;[\mathrm{mol}/\mathrm{m}^3]$",'Interpreter','latex')

	plots=get(gca, 'Children');
	l = legend(plots([20,19,18,17,4,16,15,14,13,3,12,11,10,9,2,8,7,6,5,1]),'NumColumns',4,'Interpreter','latex','FontSize',8);
	l.Layout.Tile = 'south';
	f.Units = "centimeters";
	f.Position = [8,1,19,15];
	uistack(plots(8),'top')
	uistack(plots(12),'top')
	uistack(plots(16),'top')
	uistack(plots(20),'top')
	uistack(plots(4),'top')
	uistack(plots(3),'top')
	uistack(plots(2),'top')
	uistack(plots(1),'top')

	figure(500002)
		print(gcf, "./Figures/CLAvarage.png",'-dpng','-r1200')
		print(gcf, "./Figures/CLAvarage.jpg",'-djpeg','-r1200')
		print(gcf, "./Figures/CLAvarage.eps",'-depsc','-r1200')

end
