function [physics, tvec, CL_vec]  = main_Static(model, l, u)
	fprintf('Starting job:'+model+'_'+string(l)+'_'+string(u)+'\n');

	%model = "Subgrid";  %Subgrid  WuLorenzis
	% l = 0.25e-3;
	% u = 0.01e-3;

	k = [1e-4,	1e-10,	0.5,	0;
	   	1e-10,	0,		0.3,	0;
	   	1e-6,	0,		0,		0;
	   	1e1,	7e5,	0,		0; 
	   	1e-8,	1e-13,	0.5,	0;
	   	1e-10,	1e-14,	0.3,	0;
	   	3e-5/(2*96485.3329),3e-5/(2*96485.3329), 0.5, -0.4]; %3e-7
	sname = model+'_'+string(l)+'_'+string(u);

	
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
	nmax = 1*(length(Files)-2)-2;
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
    	% mesh properties
    	mesh_in.type = "Square";
    	mesh_in.Nx    = 25;
    	mesh_in.Ny    = 250;
    	mesh_in.Lx    = 5e-3;
    	mesh_in.Ly    = 50e-3;
    	mesh_in.nfrac = 4;
    	mesh_in.ipcount1D = 4;
    	mesh_in.zeroWeight = false;
	
        dofs_in.dofs = {"dx","dy","phi","CL","Epot", "theta","H", "OH","Na", "Cl","Fe","FeOH"};
        dofs_in.Step = [2,   2,    1,    3,   3,      3,      3,   3,   3,    3,   3,    3];

    	%physics models
    	physics_in{1}.type = "PhaseFieldDamage";
    	physics_in{1}.Egroup = "Internal";
    	physics_in{1}.young = 200e9;
    	physics_in{1}.poisson = 0.3;
		physics_in{1}.kmin = 1e-10;
		sc = 3e+08*100;   %l=0.5e-3; Gc = 1e3;
		%sc = 9/16*sqrt(physics_in{1}.young*physics_in{1}.Gc/(3*physics_in{1}.l))
		physics_in{1}.l=l;
		physics_in{1}.Gc=3*16^2/9^2*physics_in{1}.l/physics_in{1}.young*sc^2;
		physics_in{1}.GDegrade = 0.9;
		physics_in{1}.NL = 1e6;
		physics_in{1}.gb = 30e3;
		sc = 9/16*sqrt(physics_in{1}.young*physics_in{1}.Gc/(3*physics_in{1}.l))
	
    	physics_in{2}.type = "HydrogenDiffusion";
		physics_in{2}.Egroup = "Internal";
    	physics_in{2}.DL = 1e-9;
		physics_in{2}.NL = 1e6;
		physics_in{2}.gb = physics_in{1}.gb;
		physics_in{2}.NT = 1e2;
		physics_in{2}.kmin = physics_in{1}.kmin;
	
    	physics_in{3}.type = "Constrainer";
    	physics_in{3}.Ngroup = "Bottom";
    	physics_in{3}.dofs = {"dy"};
    	physics_in{3}.conVal = [0];
	
		physics_in{4}.type = "Constrainer";
    	physics_in{4}.Ngroup = "LeftBottom";
    	physics_in{4}.dofs = {"dx"};
    	physics_in{4}.conVal = [0];
	
    	physics_in{5}.type = "Constrainer";
    	physics_in{5}.Ngroup = "Top";
    	physics_in{5}.dofs = {"dy"};
    	physics_in{5}.conVal = [u];

% 		physics_in{7}.type = "DummyPhaseField";
% 		physics_in{7}.l=physics_in{1}.l;
% 		physics_in{7}.Egroup = "Internal";
	
		physics_in{6}.type = "PhaseFieldElectrolyte";
    	physics_in{6}.Egroup = "Internal";
    	physics_in{6}.D = [9.3; 5.3; 1.3; 2; 1.4; 1]*1e-9;  %H OH Na CL Fe FeOH
		physics_in{6}.z = [1; -1; 1; -1; 2; 1];
		physics_in{6}.pH0 = 5;
		physics_in{6}.NaCl = 0.6e3;
		physics_in{6}.Lumped = [true; true]; %water, metal
		physics_in{6}.k = [1e6; 1e-1; 1e-3; 1e-3]; %water, Fe, Fe', FeOH
		physics_in{6}.NAds = 1e-3;
		physics_in{6}.ksurf = k;
		physics_in{6}.NL = physics_in{1}.NL;
		physics_in{6}.Em = 0;
		physics_in{6}.Lumpedsurf = [1 1 1 1 1 1 1];
		physics_in{6}.h0 = 1e-12;
		physics_in{6}.Flowtype = model;
		physics_in{6}.l = physics_in{1}.l;
		physics_in{6}.phi_trigger = -0.05;

		physics_in{7}.type = "Constrainer";
    	physics_in{7}.Ngroup = "LeftTop";
    	physics_in{7}.dofs = {"dx"};
    	physics_in{7}.conVal = [0];
	
	
    	%% solver inputs
    	solver_in.maxIt = 100;
    	solver_in.Conv = 1e-6;
    	solver_in.tiny = 1e-4;
    	solver_in.linesearch = false;
    	solver_in.linesearchLims = [0.1 1];
		solver_in.OuterLoops = 10;
	
    	%% initialization
    	mesh = Mesh(mesh_in);
    	%mesh.plot(true, true, true);
    	mesh.check();
	
    	physics = Physics(mesh, physics_in, dofs_in);
	
	
    	dt = 30;  %30
    	physics.time = 0;
	
    	n_max = 1000*24*360;
    	solver = Solver(physics, solver_in);
    	tvec = 0;
		CL_vec = 0;
		Cmax_vec = 0;
		LFrac = 0;
	
		tmax = 60*60*200;
	
    	startstep = 1;
	else
    	filename = savefolder+string(restart_num);
    	load(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax","LFrac");
    	startstep = restart_num+1;
		solver.linesearch = false;
		solver.linesearchLims = [0.05 1.0];
	end
	
	if plotonly
    	plotres(physics)
	else
    	for tstep = startstep:n_max
        	disp("Step: "+string(tstep));
			disp("Time: "+string(physics.time));
%  			if tstep==1
%  				physics.dt = 1e-3;
%  			else
				physics.dt = dt*1.05^(tstep-1);
			%end
			disp("dTime: "+string(physics.dt));
        	
        	solver.Solve();
        	
        	physics.time = physics.time+physics.dt;
        	tvec(end+1) = tvec(end)+physics.dt;
			CL_vec(end+1) = physics.models{2}.CL_int./mesh.Area(1);
			Cmax_vec(end+1) = physics.models{2}.CL_max;
			LFrac(end+1) = physics.models{1}.LFrac;
    	
        	%close all
        	if mod(tstep, 1) == 0
				plotres(physics, tvec, CL_vec, LFrac);
			end
        	if mod(tstep, 1) == 0
            	filename = savefolder+string(tstep);
            	save(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax","LFrac");
			end
	
			if (physics.time>tmax)
				break
			end
		end
		filename = savefolder+"end";
    	save(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax","LFrac");
	end
	
	toc(tmr)
end

function plotres(physics, tvec, CL_vec, LFrac)
    figure(42)
	clf 
    subplot(2,3,1)
        physics.PlotNodal("phi",0, "Internal");
        title("\phi")
		colorbar
	subplot(2,3,2)
        physics.PlotIP("sh","Internal");
        title("s_H")
		colorbar
	subplot(2,3,3)
        physics.PlotNodal("CL",0, "Internal");
        title("C_L")
		colorbar
	subplot(2,3,4)
		plot(tvec/3600, CL_vec)
		xlabel('t [hours]')
		ylabel('avarage C_L [mol/m^3]')
    subplot(2,3,5)
        physics.PlotNodal("CL",1000, "Internal");
        title("C_L")
		colorbar
	subplot(2,3,6)
		plot(tvec/3600, LFrac*1000)
		xlabel('t [hours]')
		ylabel('L_{frac} [mm]')

	figure(43)
		clf
		physics.models{6}.plotFields(physics);

     drawnow();
end

