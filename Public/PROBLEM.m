classdef PROBLEM < handle
    properties(SetAccess = private)
        Global; % The current GLOBAL object
    end
    
    properties(SetAccess = ?PROBLEM)
        Parameter;
    end
    methods(Access = protected)
        %% Constructor
        function obj = PROBLEM()
            obj.Global = GLOBAL.GetObj();
        end
    end
    
    methods
        
    end
end