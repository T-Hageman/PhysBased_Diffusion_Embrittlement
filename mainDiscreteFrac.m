function [physics, tvec, CL_vec]  = mainDiscreteFrac(irun)

	Uextern = 10^(-8+irun);   % x5

	Em = 0;

	fprintf('Starting job:'+string(Em)+'\n');

	k = [1e-4,	1e-10,	0.5,	0;
	   	1e-10,	0,		0.3,	0;
	   	1e-6,	0,		0,		0;
	   	1e1,	7e5,	0,		0; 
	   	1e-8,	1e-13,	0.5,	0;
	   	1e-10,	1e-14,	0.3,	0;
	   	3e-5/(2*96485.3329),3e-5/(2*96485.3329), 0.5, -0.4]; 
	sname = "Discrete_"+string(Uextern);
	
	maxNumCompThreads(8);
	savefolder = "./Results_NoProp/"+sname;
	mkdir(savefolder);
	savefolder=savefolder+"/";
	
	addpath(genpath('./Models'))
	addpath(genpath('./Shapes'))
	
	delete(gcp('nocreate'))
	parpool('threads')
	
	tmr = tic;
	
	%% input properties
	
	Files=dir(savefolder);
	nmax = 10*(length(Files)-2);
	restartfrom = nmax;

 	if restartfrom>0
 		restart = true;
 		restart_num = restartfrom;
 	else
		restart = false;
		restart_num = 0;
	end

	plotonly = false;
	if restart == false
        dofs_in.dofs = {"dx","dy","CL","Epot", "theta","H", "OH","Na", "Cl","Fe","FeOH"};
        dofs_in.Step = [2,   2,     1,   1,      1,      1,   1,   1,    1,   1,    1];

    	% mesh properties
    	mesh_in.type = "Square_WithDisc";
    	mesh_in.Nx    = 25;
    	mesh_in.Ny    = 250;
    	mesh_in.Lx    = 5e-3;
    	mesh_in.Ly    = 50e-3;
		mesh_in.Uextern = Uextern;
    	mesh_in.ipcount1D = 3;
    	mesh_in.zeroWeight = false;
	
    	%physics models
    	physics_in{1}.type = "LinearElastic";
    	physics_in{1}.Egroup = "Metal";
    	physics_in{1}.young = 200e9;
    	physics_in{1}.poisson = 0.3;
	
    	physics_in{2}.type = "HydrogenDiffusionUnDamaged";
		physics_in{2}.Egroup = "Metal";
    	physics_in{2}.DL = 1e-9;
		physics_in{2}.NL = 1e6;
		physics_in{2}.gb = 30e3;
		physics_in{2}.NT = 1e2;
	
    	physics_in{3}.type = "Constrainer";
    	physics_in{3}.Ngroup = "M_Bottom";
    	physics_in{3}.dofs = {"dx"};
    	physics_in{3}.conVal = [0];
	
    	physics_in{4}.type = "Constrainer";
    	physics_in{4}.Ngroup = "M_Bottom";
    	physics_in{4}.dofs = {"dy"};
    	physics_in{4}.conVal = [0];
	
		physics_in{5}.type = "Constrainer";
    	physics_in{5}.Ngroup = "M_Top";
    	physics_in{5}.dofs = {"dx"};
    	physics_in{5}.conVal = [0];
	
    	physics_in{6}.type = "Constrainer";
    	physics_in{6}.Ngroup = "M_Top";
    	physics_in{6}.dofs = {"dy"};
    	physics_in{6}.conVal = [Uextern];
	
		physics_in{7}.type = "Electrolyte";
    	physics_in{7}.Egroup = "Electrolyte";
    	physics_in{7}.D = [9.3; 5.3; 1.3; 2; 1.4; 1]*1e-9;  %H OH Na CL Fe FeOH
		physics_in{7}.z = [1; -1; 1; -1; 2; 1];
		physics_in{7}.pH0 = 5;
		physics_in{7}.NaCl = 0.6e3;
		physics_in{7}.Lumped = [true; true]; %water, metal
		physics_in{7}.k = [1e6; 1e-1; 1e-3; 1e-3]; %water, Fe, Fe', FeOH

		physics_in{8}.type = "ElectrolyteInterface";
    	physics_in{8}.Egroup = "Interface";
		physics_in{8}.NAds = 1e-3;
		physics_in{8}.k = k;
		physics_in{8}.NL = 1e6;
		physics_in{8}.Em = Em;
		physics_in{8}.Lumped = [1 1 1 1 1 1 1];

		physics_in{9}.type = "ElectrolyteInterface";
    	physics_in{9}.Egroup = "Interface2";
		physics_in{9}.NAds = 1e-3;
		physics_in{9}.k = k;
		physics_in{9}.NL = 1e6;
		physics_in{9}.Em = Em;
		physics_in{9}.Lumped = [1 1 1 1 1 1 1];
	
		initH = 1000*10^(-physics_in{7}.pH0);
		initOH = 1000*10^(-14+physics_in{7}.pH0);
		initCl = physics_in{7}.NaCl;
		initNa = initCl-initH+initOH;
	
		physics_in{10}.type = "Constrainer";
    	physics_in{10}.Ngroup = "E_Left";
    	physics_in{10}.dofs = {"Epot";"H";"OH";"Na";"Cl";"Fe";"FeOH"};
    	physics_in{10}.conVal = [0; initH; initOH; initNa; initCl; 0; 0];
	
	
	
	
    	%% solver inputs
    	solver_in.maxIt = 100;
    	solver_in.Conv = 1e-9;
    	solver_in.tiny = 1e-9;
    	solver_in.linesearch = true;
    	solver_in.linesearchLims = [0.1 1];
		solver_in.OuterLoops = 1;
	
    	%% initialization
    	mesh = Mesh(mesh_in);
    	mesh.plot(true, true, true);
    	mesh.check();
	
    	physics = Physics(mesh, physics_in, dofs_in);
	
	
    	dt = 300;
    	physics.time = 0;
	
    	n_max = 1000*24*360;
    	solver = Solver(physics, solver_in);
    	tvec = 0;
		CL_vec = 0;
		Cmax_vec = 0;
	
		tmax = 60*60*200;
	
    	startstep = 1;
	else
    	filename = savefolder+string(restart_num);
    	load(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax");
    	startstep = restart_num+1;
		solver.linesearch = false;
		solver.linesearchLims = [0.1 1];
	end
	
	if plotonly
    	plotres(physics)
	else
    	for tstep = startstep:n_max
        	disp("Step: "+string(tstep));
			disp("Time: "+string(physics.time));
% 			if tstep==1
% 				physics.dt = 1e0;
% 			else
				physics.dt = dt*1.05^(tstep-1);
%			end
			disp("dTime: "+string(physics.dt));
        	
        	solver.Solve();
        	
        	physics.time = physics.time+physics.dt;
        	tvec(end+1) = tvec(end)+physics.dt;
			CL_vec(end+1) = physics.models{2}.CL_int./mesh.Area(1);
			Cmax_vec(end+1) = physics.models{2}.CL_max;
    	
        	%close all
        	if mod(tstep, 1) == 0
            	plotres(physics, tvec, CL_vec);
			end
        	if mod(tstep, 10) == 0
            	filename = savefolder+string(tstep);
            	save(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax");
			end
	
			if (physics.time>tmax)
				break
			end
		end
		filename = savefolder+"end";
    	save(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax");
	end
	
	toc(tmr)
end

function plotres(physics, tvec, CL_vec)
    figure(42)
	clf 
    subplot(2,3,1)
        physics.PlotNodal("H",-1, "Electrolyte");
        title("H^+")
		colorbar
	subplot(2,3,2)
        physics.PlotNodal("OH",-1, "Electrolyte");
        title("OH^-")
		colorbar
	subplot(2,3,3)
        physics.PlotNodal("CL",-1, "Metal");
        title("C_L")
		colorbar
	subplot(2,3,4)
		plot(tvec/3600, CL_vec)
		xlabel('t [hours]')
		ylabel('avarage C_L [mol/m^3]')
% 	subplot(2,3,5)
% 		physics.PlotNodal("Theta",-1, "Interface");
% 		title('\theta')
% 		colorbar
	subplot(2,3,6)
        physics.PlotNodal("Fe",-1, "Electrolyte");
        title("Fe^{2+}")
		colorbar

% 	figure(44)
% 		clf
% 		physics.models{9}.plotReactions(physics);
% 
% 	figure(43)
% 		clf
% 		physics.models{7}.plotFields(physics);

     drawnow();
end

