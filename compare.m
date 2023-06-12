close all; clear all; clc

	addpath(genpath('./Models'))
	addpath(genpath('./Shapes'))

	cspec = {"r-","g-","b-","m-","r-.","g-.","b-.","m-."};

FileSet = 1;
switch FileSet
	case 1
		files = {"./Results/Subgrid_l_0.000125/";
		 		"./Results/Subgrid_l_0.00025/";
		 		"./Results/Subgrid_l_0.000375/";
		 		"./Results/Subgrid_l_0.0005/";
 		 		"./Results/WuLorenzis_l_0.000125/";
		 		"./Results//WuLorenzis_l_0.00025/";
 		 		"./Results//WuLorenzis_l_0.000375/";
 		 		"./Results//WuLorenzis_l_0.0005/";
		 		};%+"300.mat"
end


names = {"physics-based, $\ell=0.125\;\mathrm{mm}$";
		 "physics-based, $\ell=0.25\;\mathrm{mm}$";
		 "physics-based, $\ell=0.375\;\mathrm{mm}$";
		 "physics-based, $\ell=0.5\;\mathrm{mm}$";
		 "distributed diffusion, $\ell=0.125\;\mathrm{mm}$";
		 "distributed diffusion, $\ell=0.25\;\mathrm{mm}$";
		 "distributed diffusion, $\ell=0.375\;\mathrm{mm}$";
		 "distributed diffusion, $\ell=0.5\;\mathrm{mm}$";};

snames = {"p125"; "p250"; "p375"; "p500"; "d125"; "d250"; "d375"; "d500"};

m = ["x","o","d","s"]; m = [m m];

for i=1:length(files)
	Fls=dir(files{i});
	nmax = 1*(length(Fls)-2)-2;

	if (nmax>300)
		file10days = files{i}+"300.mat";
	else
		file10days = files{i}+string(nmax)+".mat";
	end

 	load(file10days, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax","LFrac");
 	ifig = [2*i-1 2*i];
 	%plotres(physics, tvec, CL_vec, LFrac, ifig)
	%f = figure(ifig(1))
	%savename = "./Figures_Dynamic/CL_"+snames{i};
	%print(gcf, savename+".png",'-dpng','-r1200')
	%print(gcf, savename+".jpg",'-djpeg','-r1200')
	%print(gcf, savename+".eps",'-depsc','-r1200')	

	fileMax = files{i}+string(nmax)+".mat";
	load(fileMax, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax","LFrac");

	figure(1001)
		hp = plot(tvec(2:end)/3600, LFrac(2:end)*1000,cspec{i}+m(i),'LineWidth',1.5);
		plotsparsemarkers(hp, [], m(i), 10 ,false);
		xlabel('$\mathrm{time} \;[\mathrm{hours}]$','Interpreter','latex')
		ylabel('$a\;[\mathrm{mm}]$','Interpreter','latex')
		hold on

	figure(1002)
		hp = plot(tvec(2:end)/3600, CL_vec(2:end),cspec{i}+m(i),'LineWidth',1.5);
		plotsparsemarkers(hp, [], m(i), 10 ,false);
		xlabel('$\mathrm{time}\;[\mathrm{hours}]$','Interpreter','latex')
		ylabel('$\overline{C}_\mathrm{L}\;[\mathrm{mol}/\mathrm{m}^3]$','Interpreter','latex')
		hold on
end

snames = {"LFrac"; "CL"}
for i=1:2
	f = figure(1000+i);
	legend(names,'Interpreter','latex','FontSize',8,'NumColumns',2,'Location','SouthOutside');
	f.Units = "centimeters";
	f.Position = [8,1,14,8];
	if i==1 
		ylim([6 14]);
	end
	xlim([0 400])
	savename = "./Figures_Dynamic/"+snames{i};
	print(gcf, savename+".png",'-dpng','-r1200')
	print(gcf, savename+".jpg",'-djpeg','-r1200')
	print(gcf, savename+".eps",'-depsc','-r1200')
end

function plotres(physics, tvec, CL_vec, LFrac, ifig)
    f = figure(ifig(1))
	clf 
%     subplot(2,3,1)
%         physics.PlotNodal("phi",0, "Internal");
%         title("\phi")
% 		colorbar
% 	subplot(2,3,2)
%         physics.PlotIP("sh","Internal");
%         title("s_H")
% 		colorbar
% 	subplot(2,3,3)
        physics.PlotNodal("CL",0, "Internal");
		cb = colorbar
		cb.Title.String = {'$C_\mathrm{L}$', '[$\mathrm{mol}/\mathrm{m}^3$]'};
		cb.Title.Interpreter='latex';
		clim([6 9])
	f.Units = "centimeters";
	f.Position = [8,1,8,6];
	axis off
	cb.Position = [0.8 0.1 0.05 0.7];
	ax = gca;
	ax.Position = [0.05 0.05 0.7 0.9];
% 	subplot(2,3,4)
% 		plot(tvec/3600, CL_vec)
% 		xlabel('t [hours]')
% 		ylabel('avarage C_L [mol/m^3]')
%     subplot(2,3,5)
%         physics.PlotNodal("CL",1000, "Internal");
%         title("C_L")
% 		colorbar
% 	subplot(2,3,6)
% 		plot(tvec/3600, LFrac*1000)
% 		xlabel('t [hours]')
% 		ylabel('L_{frac} [mm]')
% 
% 	figure(ifig(2))
% 		clf
% 		physics.models{6}.plotFields(physics);

     drawnow();
end