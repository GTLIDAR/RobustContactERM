classdef TwiceSmoothStepTerrain < SmoothStepTerrain
    
    properties (Constant)
       spline_order = 5; 
    end
    
    methods (Access = protected)
        function obj = solveCoeff(obj)
            
            % Calculate the coefficients of the step transition polynomial
            x0 = obj.step_location;
            x1 = obj.step_location + obj.step_width;
            y0 = obj.terrain_height;
            y1 = obj.step_height;
            
            [p0, D] = obj.polynomialize(x0);
            p1 = obj.polynomialize(x1);
            A = [p0, p1, D*p0, D*p1, D^2*p0, D^2*p1]';
            b = [y0; y1; 0; 0; 0; 0];
            % Solve for the polynomial coefficients.
            obj.spline_coeff = A\b;
        end
    end
end