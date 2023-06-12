clear all
close all
clc

surfplots = false;

addpath(genpath('./Models'))
addpath(genpath('./Shapes'))
	
cspec = {"r","g","b","m"};
lspec = {"-.","--","-",":"};

i=0;
for modeli=1:2   
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

%%Combined and individual plots
figure(500000)
figure(500001)
figure(500002)
tiledlayout(1,2);
for i=1:length(simsets)
	sname = simsets{i}.model+'_'+string(simsets{i}.l)+'_'+string(simsets{i}.u);
	savefolder = "./Results_NoProp/"+sname;
	filename = savefolder+"/end";
	dispname = "$\ell="+string(simsets{i}.l*1000)+",\; U_{\mathrm{ext}}="+string(simsets{i}.u*1000)+"$";
	load(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax","LFrac");
	if (simsets{i}.model == "Subgrid")
		figure(500002)
		nexttile(1);
	else
		figure(500002)
		nexttile(2);
	end
	plot(tvec/3600,CL_vec,simsets{i}.linespec,'DisplayName',dispname,'LineWidth',2);
	hold on

	if (simsets{i}.model == "Subgrid")
		figure(500000)
	else
		figure(500001)
	end	
	plot(tvec/3600,CL_vec,simsets{i}.linespec,'DisplayName',dispname,'LineWidth',2);
	hold on

	if surfplots
		f1 = figure(i);
			physics.PlotNodal("phi",0, "Internal");
			caxis([0 1]);
			cb = colorbar;
			cb.Title.String = {'$\phi$', '[$-$]'};
			cb.Title.Interpreter='latex';
			axis off
			saveFigNow(f1, "./Figures_Loose/PhaseField_"+sname, 3, true, true, cb)
		f2 = figure(100+i);
			physics.PlotNodal("CL",0, "Internal");
			caxis([0 5]);
			cb = colorbar;
			cb.Title.String = {'$C_\mathrm{L}$', '[$\mathrm{mol}/\mathrm{m}^3$]'};
			cb.Title.Interpreter='latex';
			axis off
			saveFigNow(f2, "./Figures_Loose/CL_"+sname, 3, true, true, cb)
		f3 = figure(200+i);
			physics.models{6}.plotFieldspH(physics);
			caxis([5 6]);
			hold on
			plot([0 5e-3 5e-3 0 0],[0 0 50e-3 50e-3 0],'k-')
			cb = colorbar;
			cb.Title.String = {'$\mathrm{pH}$', '[$-$]'};
			cb.Title.Interpreter='latex';
			axis off
			saveFigNow(f3, "./Figures_Loose/pH_"+sname, 3, true, true, cb)
		
		close(f1)
		close(f2)
		close(f3)
	end
end

mrkrs = {"kx-","kd-","ko-","ks-","k-"}
figure(500003)
for j=1:4
	filename = "./Results_NoProp/"+"Discrete_"+string(10^(-8+j))+"/end";
	load(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax","LFrac");
	figure(500002)
	nexttile(1);
	plot(tvec/3600,CL_vec,mrkrs{j},'DisplayName','Discrete $U_{\mathrm{ext}}='+string(10^(-8+j+3))+"$");
	nexttile(2);
	plot(tvec/3600,CL_vec,mrkrs{j},'DisplayName','Discrete $U_{\mathrm{ext}}='+string(10^(-8+j+3))+"$");
	figure(500003)
	plot(tvec/3600,CL_vec,mrkrs{j},'DisplayName','Discrete $U_{\mathrm{ext}}='+string(10^(-8+j+3))+"$");
	hold on
	figure(500001)
	plot(tvec/3600,CL_vec,mrkrs{j},'DisplayName','Discrete $U_{\mathrm{ext}}='+string(10^(-8+j+3))+"$");
	figure(500000)
	plot(tvec/3600,CL_vec,mrkrs{j},'DisplayName','Discrete $U_{\mathrm{ext}}='+string(10^(-8+j+3))+"$");
end

for i=1:4
f = figure(500000+i-1);
	if(i==3)
		nexttile(1)
	end
	ax1 = gca;
	xlim([0 200]);
	ylim([0 1.5]);
	xlabel("$\mathrm{time}\;[hours]$",'Interpreter','latex');
	ylabel("$\overline{C}_\mathrm{L}\;[\mathrm{mol}/\mathrm{m}^3]$",'Interpreter','latex')
	if (i==3)
		nexttile(2)
		ax2 = gca;
		xlim([0 200]);
		ylim([0 1.5]);
		xlabel("$\mathrm{time}\;[hours]$",'Interpreter','latex');
		ylabel("$\overline{C}_L\;[\mathrm{mol}/\mathrm{m}^3]$",'Interpreter','latex')
	end

	if (i<3.5)
	plots=get(gca, 'Children');
	l = legend(plots([20,19,18,17,4,16,15,14,13,3,12,11,10,9,2,8,7,6,5,1]),'NumColumns',4,'Interpreter','latex','FontSize',8);
	if (i==3)
	l.Layout.Tile = 'south';
	else
	l.Location = 'southoutside';
	end
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
	else
	l = legend('NumColumns',4,'Interpreter','latex','FontSize',8);
	l.Location = 'southoutside';
	end

	saveFigNow(f, "./Figures_Loose/CLAvarage"+string(i), 19, true, false, 0)
end


function saveFigNow(fg, sname, HFig, SaveDouble, hasColorbar, cb)
	fprintf(sname+"  ")

	fg.Units = 'centimeters';
	if (SaveDouble)
		fg.Position = [2 2 19 HFig];
	else
		fg.Position = [2 2 8 HFig];
	end
	if (hasColorbar)
		cb.Position(1) = 0.9;
		cb.Position(4) = 0.6;
		ax = gca;
		ax.Position(1) = 0.05;
		ax.Position(3) = 0.8;
	end

	drawnow();
	print(fg, sname+".png",'-dpng','-r1200'); fprintf(".png  ")
	print(fg, sname+".jpg",'-djpeg','-r1200'); fprintf(".jpg  ")
	print(fg, sname+".eps",'-depsc','-r1200'); fprintf(".eps  ")
	print(fg, sname+".emf",'-dmeta','-r1200'); fprintf(".emf\n")
end




