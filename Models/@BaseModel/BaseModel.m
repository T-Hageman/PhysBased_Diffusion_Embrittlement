classdef BaseModel < handle
    %CONSTRAINER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function obj = BaseModel()

        end
        
        function Commit(obj, physics, commit_type)
            
        end
        
        function Irr = Irreversibles(obj, physics)
            Irr = false;
        end
        
        function [hasInfo, provided] = Provide_Info(obj, physics, var, elems, loc)
           hasInfo = false;
           provided = [];
        end
        
        function getKf(obj, physics, stp)

		end

		function OncePerStep(obj, physics, stp)

		end
    end
end

