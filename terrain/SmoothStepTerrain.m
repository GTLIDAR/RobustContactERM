classdef SmoothStepTerrain < Terrain2D
    %SMOOTHSTEPTERRAIN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        friction_coeff(1,1) double {mustBeFinite, mustBeReal, mustBeNonnegative} = 0.5;
        terrain_height(1,1) double {mustBeFinite, mustBeReal} = 0.0;
        step_location(1,1) double {mustBeReal} = 0.5;
        step_height(1,1) double {mustBeReal} = 0.5;
        step_width(1,1) double {mustBeReal, mustBePositive} = 0.05;
    end
    properties (Hidden)
        % spline_coeff: The coefficients of the polynomial spline
        % connecting the lower and upper terrains.
        spline_coeff = [];
    end
    properties (Abstract, Constant)
       spline_order; 
    end
    %% Abstract Methods: These methods define the smoothness of the terrain
    methods (Abstract, Access = protected)
       obj = solveCoeff(obj);
       % solveCoeff: solve for the polynomial coefficients of the step
       % transition
    end
    %% Methods common to all implementations
    methods
        function obj = SmoothStepTerrain()
            %% SmoothStepTerrain: Step Terrains with Smooth transtiions
            %
            %  SmoothStepTerrain is an abstract class for implementing a
            %  step terrain with a smooth transition. Concrete classes can
            %  implement SmoothStepTerrain by specifying the transition as
            %  a polynomial spline.
            %
            %   See also: OnceSmoothStepTerrain, TwiceSmoothStepTerrain
           
            % Solve for the polynomial spline coefficients
            obj  = obj.solveCoeff();
        end
        function x = nearest_opt(obj, xA)
           x = fmincon(@(x) sum((xA-x).^2), xA, [], [],[], [], [], [], @(x)obj.terrainWrapper(x)); 
        end
        function y = splineEval(obj, x)
            y = obj.spline_coeff' * obj.polynomialize(x);
        end
        function [c, ceq] = splineCstr(obj,x)
           y = obj.spline_coeff' * obj.polynomialize(x(1));
           ceq = x(2) - y;
           c = 0;
        end
        function x = spline_nearest(obj, xA)
           x = fmincon(@(x) (x(1) - xA(1))^2 +(x(2)-xA(2))^2, xA, [], [], [], [], [], [], @(x)obj.splineCstr(x)); 
        end
        function [x, dx] = nearest(obj,xA)
            %% Nearest: Returns the nearest point on the terrain
            
            dim = numel(xA);
            % Candidate points on the polynomial spline
            % Solve the polynomial rootfinding problem
            coeff = obj.spline_coeff';
            [~,D] = obj.polynomialize(xA(1));
            dcoeff = coeff*D;
            
            coeff(end) = coeff(end) - xA(2);
            
            pol = conv(coeff,dcoeff);
            pol(end-1) = pol(end-1)+1;
            pol(end) = pol(end) - xA(1);
            % Canditate points using rootfinding
            x_cand = roots(pol);
            % Filter out complex roots
            x_imag = imag(x_cand);
            x_cand(abs(x_imag) > 1e-12) = [];
            x_cand = real(x_cand);
            
            % Filter out the roots which are out of range
            x_cand(or(x_cand < obj.step_location, x_cand > obj.step_location + obj.step_width)) = [];
            if ~isempty(x_cand)
                y_cand = obj.polynomialize(x_cand')' * obj.spline_coeff;
            else
                y_cand = [];
            end
            % The nearest point on the terrain
            if xA(1) < obj.step_location
                x_cand = [xA(1); x_cand; obj.step_location + obj.step_width];
            elseif xA(1) > obj.step_location + obj.step_width
                x_cand = [obj.step_location; x_cand; xA(1)];
            else
                x_cand = [obj.step_location; x_cand; obj.step_location + obj.step_width];
            end
            y_cand = [obj.terrain_height; y_cand; obj.step_height];
                        
            dist = (xA(1) - x_cand).^2 + (xA(2) - y_cand).^2;
            [~,idx] = min(dist);
            
            x = [x_cand(idx); y_cand(idx)];
            
            if nargout > 1
                [~, dg, Hg] = obj.terrainEval(x);
                % The Lagrange Multiplier
                lambda = (xA(2) - x(2));
                % The derivative of the nearest point
                H = [eye(dim) + lambda*Hg, dg; dg', 0];
                b = [eye(dim); zeros(1,dim)];
                % Derivative of the nearest point
                dr = H\b;
                % Separete the derivative of the lagrange multiplier from the
                % derivative of the point
                %dlambda = dr(dim+1:end,:);
                dx = dr(1:dim,:);
            end
        end
        function [c, ceq] = terrainWrapper(obj,x)
           c = 0;
           ceq = obj.terrainEval(x);
        end
        function [g, dg, Hg] = terrainEval(obj, x)
            
            if x(1) < obj.step_location
                % Lower Step
                g = x(2) - obj.terrain_height;
                dg = [0;1];
                Hg = zeros(2);
                
            elseif x(1) > obj.step_location + obj.step_width
                % Upper step
                g = x(2) - obj.step_height;
                dg = [0;1];
                Hg = zeros(2);
            else
                % Spline terrain
                [p, D] = obj.polynomialize(x(1));
                g = x(2) - obj.spline_coeff'*p;
                dg = [-obj.spline_coeff'*D*p; 1];
                Hg = zeros(2);
                Hg(1,1) = -obj.spline_coeff'*D^2*p;
            end
        end
        function [N,T,dN,dT] = basis(obj, x)
            %% Basis
            
            
            if x(1) < obj.step_location || x(1) > obj.step_location + obj.step_width
                % Normal and Tangent vectors for the flat part of the terrain
                N = [0,1]';
                dN = zeros(2,2);
                T = [1,0]';
                dT = zeros(2,2);
            else
                % Normal and tangent vectors for the curved part of the
                % terrain
                % Calculate the derivatives of the polynomial part:
                [p,D] = obj.polynomialize(x(1));
                % Evaluate the polynomial and calculate the derivatives
                dp = obj.spline_coeff'*D*p;
                ddp = obj.spline_coeff'*D^2*p;
                % Now calculate the normal and tangential vectors
                N = [-dp,1]';
                T = [1, dp]';
                % Derivatives of the normal and tangent vectors
                dN = diag([-ddp, 0]);
                dT = zeros(2,2);
                dT(2,1) = ddp;
                % Normalize the Normal Vector
                sN = norm(N);
                N = N./sN;
                % Normalized Derivative
                dN = (eye(numel(N)) - N*N')*dN./sN;
                
                % Normalize the Tangent Vector
                sT = norm(T);
                T = T./sT;
                % Normalized Derivative
                dT = (eye(numel(T)) - T*T')*dT./sT;
            end
        end
        function [x,y] = draw(obj,xlim, N)
            %% Draw
            
            x = linspace(xlim(1), xlim(2), N);
            y = zeros(1,N);
            % Set the level parts of the terrain
            y(x <= obj.step_location) = obj.terrain_height;
            y(x >= obj.step_location + obj.step_width) = obj.step_height;
            % Create the polynomial parts of the terrain
            smooth_idx = and(x> obj.step_location, x<obj.step_location + obj.step_width);
            if any(smooth_idx)
                p = obj.polynomialize(x(smooth_idx));
                y(smooth_idx) = obj.spline_coeff'*p;
            end
            % If there are no outputs, plot the result.
            if nargout == 0
               figure();
               plot(x,y,'k-','LineWidth',1.5);
            end
        end
    end
    methods
        function obj = set.terrain_height(obj, val)
            obj.terrain_height = val;
            obj = obj.solveCoeff();
        end
        function obj = set.step_location(obj, val)
            obj.step_location = val;
            obj = obj.solveCoeff();
        end
        function obj = set.step_height(obj, val)
            obj.step_height = val;
            obj = obj.solveCoeff();
        end
        function obj = set.step_width(obj, val)
            obj.step_width = val;
            obj = obj.solveCoeff();
        end
    end
    methods 
        function  [p, D] = polynomialize(obj, x)
            % polynomialize: return the first N powers of the input scalar x,
            % along with a derivative operator matrix D such that D^n*p is the
            % nth derivative of the powers of x with respect to x
            v = (obj.spline_order:-1:0)';
            p = x.^v;
            D = diag(v(1:end-1),1);
        end
    end
    methods (Static)
        function [g, dg] = distance(p,r)
           g = norm(p-r);
           dg = 2*(p-r)';
        end
    end
end

