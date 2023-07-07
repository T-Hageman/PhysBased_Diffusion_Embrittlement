classdef Solver < handle
    %SOLVER Holds data pertaining to the solving of a nonlinear system of
	%equations. Required input data:
	%   solver_in.maxIt = 100;	% Maximum number of Newton-Raphson iterations per time step
    %	solver_in.Conv = 1e-9;	% Convergence criterion (relative energy based)
    %	solver_in.tiny = 1e-9;	% Convergence criterion (absolute energy based)
    %	solver_in.linesearch = true; % Flag to allow for a linear line-search
    %	solver_in.linesearchLims = [0.1 1]; % Bounding box for line-search algorith
    
    properties
        physics		%pointer to the physics object defining the problem
		NSteps		%number of staggered steps
        
        maxIt		%maximum number of iterations
        Conv		%relative convergence criterium
        tiny		%factor below which the absolute convergence criterion is fulfilled
        linesearch	%flag to enable linesearch
        linesearchLims	%limits in which the line-search algorithm operates

		convergence_log	%logged convergence data

		OuterLoops	%max number of outer loops
    end
    
    methods
        Solve(obj); %external file, staggered NR solving procedure
        
        function obj = Solver(physics, inputs)
            obj.physics = physics;
            
            obj.maxIt = inputs.maxIt;
            obj.Conv = inputs.Conv;
            obj.tiny = inputs.tiny;
            obj.linesearch = inputs.linesearch;
			if (obj.linesearch)
                obj.linesearchLims = inputs.linesearchLims;
			end

			obj.NSteps = obj.physics.dofSpace.NSteps;
			obj.convergence_log = zeros(0,obj.NSteps,0);
			obj.OuterLoops = inputs.OuterLoops;
        end
        
	
    end
end

