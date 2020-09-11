classdef optimizationHandler < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        plant;
        options;
        solverSettings;
        stateConstraints = {};
        schedule = {};
        costweights;
        optimizer = "ROBUST";
        customCosts = {};
    end
    
    methods
        function obj = optimizationHandler(model)
            % Add the plant to the handler
            obj.plant = model;
            % Get the default settings for optimization
            obj.options.duration = [1,1];
            obj.options.meshsize = 101;
        end
        
        function addSchedule(obj, varname, values)
            %% addSchedule: Adds variables to the optimization schedule
            %
            %   when a schedule is present, optimizationHandler runs
            %   multiple optimization problems at once, each time changing
            %   a value in optimization according to the current value in
            %   the schedule. Each subsequent optimization starts from the
            %   previous solution
            
            % Check that the schedule values are valid
            
            
            % Create a new schedule using the input values
            newschedule = cell(numel(values),2);
            for n = 1:numel(values)
               newschedule(n,:) = {varname, values(n)}; 
            end
            % Append the new schedule to the existing schedule
            obj.schedule = [obj.schedule; newschedule];
        end  
        function setOptions(obj, varargin)
            
            [obj.options, unused] = setOptimizerOptions(obj.options, varargin{:});
            obj.solverOptions = setSolverOptions(obj.options,unused{:});
        end
        function obj = runOptimization(obj)
            
            
            
            
            for n = 1:length(obj.schedule)
                obj.setOptions(obj.schedule{n,1},obj.schedule{n,2});
                soln = optimizePlant(obj.plant, guess, obj.options, obj.solverOptions);
                
            end
        end
        function saveHandler(obj)
            
        end
        function addCustomCost(obj, costHandle, costName)
           obj.customCosts = [obj.customCosts; costHandle]; 
        end
    end
    methods
        function report()
            
        end
    end    
end

