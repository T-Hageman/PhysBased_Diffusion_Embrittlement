function [physics, tvec, CL_vec]  = main(irun)
	%% MATLAB code which replicates the results from XXXXXX, section 4.2 (Figures 9-11)				

	%Select which results to replicate, either using the physics-based
	%diffusion model (referred to within code as "subgrid") or the
	%distributed diffusion model (using "WuLorenzis"). Also selects the
	%phase-field length scale "l" that is being used
	if irun<4.5
		model = "Subgrid";
		l = 0.125e-3*irun;
	else
		model = "WuLorenzis";
		l = 0.125e-3*(irun-4);
	end
	sname = model+"_l_"+string(l);
	
	hpc = false; %if running on hpc, use more cpu cores but surpress figure outputs
	if (hpc)
		maxNumCompThreads(32);
	else
		maxNumCompThreads(8);
	end

	%create folder into which to save results
	savefolder = "./Results/"+sname;
	mkdir(savefolder);
	savefolder=savefolder+"/";
	
	%add required subfolders to path
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

	if restart == false
		%Start simulations from the beginning

    	% mesh properties
    	mesh_in.type = "Square";
    	mesh_in.Nx    = 100;     %number of elements in horizontal direction
    	mesh_in.Ny    = 200;     %number of elements in vertical direction
    	mesh_in.Lx    = 10e-3;   %domain length [m]
    	mesh_in.Ly    = 10e-3;   %domain height [m]
    	mesh_in.ipcount1D = 4;	 %number of gauss integration points per direction 
    	mesh_in.zeroWeight = false; %add additional zero-weight integration points on element boundary
	
		%Define used degrees of freedom, and in which staggered solution
		%step they are solved for
        dofs_in.dofs = {"dx","dy","phi","CL","Epot", "theta","H", "OH","Na", "Cl","Fe","FeOH"};
        dofs_in.Step = [2,   2,    1,    3,   3,      3,      3,   3,   3,    3,   3,    3];

    	%% physics models

		%Momentum balance and phase-field evolution
    	physics_in{1}.type = "PhaseFieldDamage";
    	physics_in{1}.Egroup = "Internal";
    	physics_in{1}.young = 200e9;	%Youngs Modulus [Pa]
    	physics_in{1}.poisson = 0.3;	%Poisson ratio [-]
		physics_in{1}.kmin = 1e-10;		%residual stiffness factor [-]
		physics_in{1}.l=l;				%Phase-field length scale [m]
		physics_in{1}.Gc=2.0e3;			%Fracture release energy [J/m^2]
		physics_in{1}.GDegrade = 0.9;	%Maximum hydrogen degradation factor
		physics_in{1}.NL = 1e6;			%Concentration of interstitial lattice sites [mol/m^3]
		physics_in{1}.gb = 30e3;		%Grain boundary binding energy [J/mol]
	
		%Interstitial lattice hydrogen diffusion model
    	physics_in{2}.type = "HydrogenDiffusion";
		physics_in{2}.Egroup = "Internal";
    	physics_in{2}.DL = 1e-9;				%Diffusivity [m/s]
		physics_in{2}.NL = physics_in{1}.NL;	%Concentration of interstitial lattice sites [mol/m^3]	
		physics_in{2}.gb = physics_in{1}.gb;	%Grain boundary binding energy [J/mol]
		physics_in{2}.NT = 1e2;					%concentration of trapping sites
		physics_in{2}.kmin = physics_in{1}.kmin;%residual stiffness factor [-]
	
		%displacement constrain at bottom boundary
    	physics_in{3}.type = "Constrainer";
    	physics_in{3}.Ngroup = "Bottom";
    	physics_in{3}.dofs = {"dy"};
    	physics_in{3}.conVal = [0];
	
		%displacement constrain at left-bottom corner
		physics_in{4}.type = "Constrainer";
    	physics_in{4}.Ngroup = "LeftBottom";
    	physics_in{4}.dofs = {"dx"};
    	physics_in{4}.conVal = [0];
	
		%displacement constrain at top boundary
    	physics_in{5}.type = "Constrainer";
    	physics_in{5}.Ngroup = "Top";
    	physics_in{5}.dofs = {"dy"};
    	physics_in{5}.conVal = [0.01e-3];	%imposed displacement [m]
	
		%Crack-contained electrolyte diffusion, electro-migration, and reactions
		physics_in{6}.type = "PhaseFieldElectrolyte";
    	physics_in{6}.Egroup = "Internal";
    	physics_in{6}.D = [9.3; 5.3; 1.3; 2; 1.4; 1]*1e-9;  %Diffusion coefficients for [H OH Na CL Fe FeOH] species respectively [m/s]
		physics_in{6}.z = [1; -1; 1; -1; 2; 1];				%ionic charges for [H OH Na CL Fe FeOH] [-]
		physics_in{6}.pH0 = 5;								%Initial and boudnary pH [-]
		physics_in{6}.NaCl = 0.6e3;							%Initial and boundary Cl- concentration [mol/m^3]
		physics_in{6}.Lumped = [true; true];				%Flag indicationg whtehr to ude lumped integration for the water autoionisation and metal ion reactions
		physics_in{6}.k = [1e6; 1e-1; 1e-3; 1e-3];			%Dummy constant for the water auto-ionisation reaction, and reaction rates for Fe, Fe', FeOH reactions
		physics_in{6}.NAds = 1e-3;							%Concentration of surface adsorption sites
		physics_in{6}.ksurf = [ 1e-4,	1e-10,	0.5,	0;	%Reaction constants for surface reactions, [k k' alpha E_eq] 
	   							1e-10,	0,		0.3,	0;
	   							1e-6,	0,		0,		0;
	   							1e1,	7e5,	0,		0; 
	   							1e-8,	1e-13,	0.5,	0;
	   							1e-10,	1e-14,	0.3,	0;
								3e-5/(2*96485.3329),3e-5/(2*96485.3329), 0.5, -0.4];
		physics_in{6}.NL = physics_in{1}.NL;				%Concentration of interstitial lattice sites [mol/m^3]
		physics_in{6}.Em = 0;								%Metal electric potential
		physics_in{6}.Lumpedsurf = [1 1 1 1 1 1 1];			%Flags to indicate the use of lumped integration for surface reactions
		physics_in{6}.h0 = 1e-12;							%small offset used to prevent ill-conditioned systems
		physics_in{6}.Flowtype = model;						% Model to use for electrolyte diffusion, either Subgrid  WuLorenzis
		physics_in{6}.l = physics_in{1}.l;					%Phase-field length scale	

		%displacement constrain at left-top corner
		physics_in{7}.type = "Constrainer";
    	physics_in{7}.Ngroup = "LeftTop";
    	physics_in{7}.dofs = {"dx"};
    	physics_in{7}.conVal = [0];
	
    	%% solver inputs
    	solver_in.maxIt = 100;		%maximum amount of iterations within the nonlinear solver
    	solver_in.Conv = 1e-6;		%Relative Energy-based convergence criterion
    	solver_in.tiny = 1e-4;		%Absolute Energy-based convergence criterion
    	solver_in.linesearch = false; %flag to indicate the use of a linear line-search
    	solver_in.linesearchLims = [0.1 1]; %limits within which the line-search is performed
		solver_in.OuterLoops = 10;	%Maximum number of staggeres dolution loops to perform
	
    	%% initialization
    	mesh = Mesh(mesh_in);
		if (hpc==false)
			mesh.plot(true, true, true);
		end
    	mesh.check();
	
    	physics = Physics(mesh, physics_in, dofs_in);
	
    	dt = 30; 
    	physics.time = 0;
	

    	solver = Solver(physics, solver_in);
    	tvec = 0;
		CL_vec = 0;
		Cmax_vec = 0;
		LFrac = 0;

	    n_max = 1000*24*360;
		tmax = 60*60*1000;
	
    	startstep = 1;
	else
		%restart from previously saved file
    	filename = savefolder+string(restart_num);
    	load(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax","LFrac");
    	startstep = restart_num+1;
	end

	%Perform simulation
	for tstep = startstep:n_max
    	disp("Step: "+string(tstep));
		disp("Time: "+string(physics.time));
		physics.dt = dt*1.05^(min(100,tstep-1));
		disp("dTime: "+string(physics.dt));
    	
		%solve for current time increment
    	solver.Solve();
    	
		%Save time-series outputs
    	physics.time = physics.time+physics.dt;
    	tvec(end+1) = tvec(end)+physics.dt;
		CL_vec(end+1) = physics.models{2}.CL_int./mesh.Area(1);
		Cmax_vec(end+1) = physics.models{2}.CL_max;
		LFrac(end+1) = physics.models{1}.LFrac;
	
		%Plot results
		if (mod(tstep, 1) && hpc==false) == 0
        	plotres(physics, tvec, CL_vec, LFrac);
		end

		%save output file for later post-processing (or restarting from)
		if mod(tstep, 1) == 0
        	filename = savefolder+string(tstep);
        	save(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax","LFrac");
		end

		%check if total time has been simulated
		if (physics.time>tmax)
			break
		end
	end

	%save results one last time
	filename = savefolder+"end";
	save(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax","LFrac");
	
	toc(tmr)
end

%Plot results
function plotres(physics, tvec, CL_vec, LFrac)
    figure(42)
	clf 

	%Phase-field
    subplot(2,3,1)
        physics.PlotNodal("phi",0, "Internal");
        title("\phi")
		colorbar

	%hydro-static stresses
	subplot(2,3,2)
        physics.PlotIP("sh","Internal");
        title("s_H")
		colorbar

	%Inerstitial lattice hydrogen concentration
	subplot(2,3,3)
        physics.PlotNodal("CL",0, "Internal");
        title("C_L")
		colorbar

	%Interstitial lattice hydrogen over time
	subplot(2,3,4)
		plot(tvec/3600, CL_vec)
		xlabel('t [hours]')
		ylabel('avarage C_L [mol/m^3]')

	%Inerstitial lattice hydrogen concentration + deformations (scaled by x1000)
    subplot(2,3,5)
        physics.PlotNodal("CL",1000, "Internal");
        title("C_L")
		colorbar

	%Fracture length over time
	subplot(2,3,6)
		plot(tvec/3600, LFrac*1000)
		xlabel('t [hours]')
		ylabel('L_{frac} [mm]')

	%Fracture contained electrolyte plots
	figure(43)
		clf
		physics.models{6}.plotFields(physics);

     drawnow();
end

