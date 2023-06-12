classdef Solver < handle
    %SOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        physics
		NSteps
        
        maxIt
        Conv
        tiny
        linesearch
        linesearchLims

		convergence_log

		OuterLoops
    end
    
    methods
        Solve(obj);
        
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

