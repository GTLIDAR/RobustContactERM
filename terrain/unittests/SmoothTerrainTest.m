classdef SmoothTerrainTest < matlab.unittest.TestCase 
        
    properties
        % Terrains
        linear = ContinuousStepTerrain();
        cubic = OnceSmoothStepTerrain();
        quintic = TwiceSmoothStepTerrain();
        % Test points
        x;
        labels = {'base','spline','step'};
        % Finite Difference parameters
        h = 1e-6;
        tol = 1e-6;
    end
    
    methods (TestMethodSetup)
        function createExamples(testCase)
            

            % Three test points:
            %   The first is closest to the bottom flat terrain
            %   The second is closest to the spline terrain
            %   The third is closest to the top flat terrain
            testCase.x = {[0.1, 0.1]'; [0.51, 0.1]'; [0.6, 0.6]'};
        end
    end
    
    methods (Test)       
        function cubicNearestTest(testCase)
            % Import the necessary packages
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.IsEqualTo
            for n = 1:length(testCase.x)
                % Test value and calcualted derivatives
                r = testCase.x{n};
                z = testCase.cubic.nearest(r);
                z_est = testCase.cubic.nearest_opt(r);
                % Check the first derivative
                testCase.verifyThat(z, IsEqualTo(z_est, 'Within', RelativeTolerance(testCase.tol) | AbsoluteTolerance(testCase.tol)),sprintf('Nearest point for cubic terrain failed on %s terrain',testCase.labels{n}));
            end
            
        end
        function quinticNearestTest(testCase)
            % Import the necessary packages
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.IsEqualTo
            for n = 1:length(testCase.x)
                % Test point and calculated derivatives
                r = testCase.x{n};
                z = testCase.quintic.nearest(r);
                z_est = testCase.quintic.nearest_opt(r);
                testCase.verifyThat(z, IsEqualTo(z_est,'Within',RelativeTolerance(testCase.tol)|AbsoluteTolerance(testCase.tol)), sprintf('Nearest point for quintic terrain failed on %s terrain',testCase.labels{n}));
            end
        end
        function cubicPolynomializeTest(testCase)
            % Import the necessary packages
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.IsEqualTo
            r = testCase.x{2};
            [p,D] = testCase.cubic.polynomialize(r(1));
            dp = D*p;
            dp_est = testCase.finiteDifference(@(z)testCase.cubic.polynomialize(z), r(1));
            testCase.verifyThat(dp, IsEqualTo(dp_est, 'Within',  RelativeTolerance(testCase.tol) | AbsoluteTolerance(testCase.tol)),sprintf('Polynomial derivatives failed on cubic terrain'));
                       
        end
        function quinticPolynomializeTest(testCase)
            % Import the necessary packages
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.IsEqualTo
            r = testCase.x{2};
            [p,D] = testCase.quintic.polynomialize(r(1));
            dp = D*p;
            dp_est = testCase.finiteDifference(@(z)testCase.quintic.polynomialize(z), r(1));
            testCase.verifyThat(dp, IsEqualTo(dp_est, 'Within',  RelativeTolerance(testCase.tol) | AbsoluteTolerance(testCase.tol)),sprintf('Polynomial derivatives failed on quintic terrain'));
            
        end
        function splineDerivativesTest(testCase)
            % Import the necessary packages
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.IsEqualTo
            r = testCase.x{2};
            r = r(1);
            
            % Cubic Spline eval
            [~,dg,Hg] = testCase.cubic_eval(r);
            dg_est = testCase.finiteDifference(@(z)testCase.cubic_eval(z), r);
            testCase.verifyThat(dg, IsEqualTo(dg_est, 'Within', RelativeTolerance(testCase.tol) | AbsoluteTolerance(testCase.tol)),sprintf('Cubic spline first derivative inaccurate'));
            hess = @(z)testCase.outputWrapper2(@(w)testCase.cubic_eval(w),z);
            Hg_est = testCase.finiteDifference(@(z)hess(z),r);
            testCase.verifyThat(Hg, IsEqualTo(Hg_est, 'Within', RelativeTolerance(testCase.tol) | AbsoluteTolerance(testCase.tol)), sprintf('Cubic spline second derivative inaccurate'));
            % Quintic Spline Eval
            [~,dg,Hg] = testCase.quintic_eval(r);
            dg_est = testCase.finiteDifference(@(z)testCase.quintic_eval(z), r);
            testCase.verifyThat(dg, IsEqualTo(dg_est, 'Within', RelativeTolerance(testCase.tol) | AbsoluteTolerance(testCase.tol)),sprintf('Quintic spline first derivative inaccurate'));
            hess = @(z)testCase.outputWrapper2(@(w)testCase.quintic_eval(w),z);
            Hg_est = testCase.finiteDifference(@(z)hess(z),r);
            testCase.verifyThat(Hg, IsEqualTo(Hg_est, 'Within', RelativeTolerance(testCase.tol) | AbsoluteTolerance(testCase.tol)), sprintf('Quintic spline second derivative inaccurate'));
        end
        function cubicNearestDerivativeTest(testCase)
            % Import the necessary packages
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.IsEqualTo
            for n = 1:length(testCase.x)
                % Test value and calcualted derivatives
                r = testCase.x{n};
                [~,dz] = testCase.cubic.nearest(r);
                % Check the first derivative
                dz_est = testCase.finiteDifference(@(z)testCase.cubic.nearest(z),r);
                testCase.verifyThat(dz, IsEqualTo(dz_est, 'Within', RelativeTolerance(testCase.tol) | AbsoluteTolerance(testCase.tol)),sprintf('Nearest derivative for cubic terrain failed on %s terrain',testCase.labels{n}));
            end
        end
        function quinticNearestDerivativeTest(testCase)
            % Import the necessary packages
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.IsEqualTo
            for n = 1:length(testCase.x)
                % Test point and calculated derivatives
                r = testCase.x{n};
                [~,dz] = testCase.quintic.nearest(r);
                % Check the first derivative
                dz_est = testCase.finiteDifference(@(z)testCase.quintic.nearest(z),r);
                testCase.verifyThat(dz, IsEqualTo(dz_est,'Within',RelativeTolerance(testCase.tol)|AbsoluteTolerance(testCase.tol)), sprintf('Nearest derivative for quintic terrain failed on %s terrain',testCase.labels{n}));
            end
        end
        function cubicBasisDerivativeTest(testCase)
            % Import the necessary packages
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.IsEqualTo
            % Create a handle to the basis function for calculating the
            % hessian
            tangentfun =  @(w)testCase.outputWrapper2(@(z)testCase.cubic.basis(z), w);
            for n = 1:length(testCase.x)
                % Test point and calculated derivatives
                r = testCase.x{n};
                [~, ~, dN, dT] = testCase.cubic.basis(r);
                
                % Check the derivative of the normal and tangent vectors
                dN_est = testCase.finiteDifference(@(z)testCase.cubic.basis(z), r);
                dT_est = testCase.finiteDifference(@(z)tangentfun(z), r);
                
                testCase.verifyThat(dN, IsEqualTo(dN_est,'Within',RelativeTolerance(testCase.tol)|AbsoluteTolerance(testCase.tol)), sprintf('Normal vector derivative for cubic terrain failed on %s terrain', testCase.labels{n}));
                testCase.verifyThat(dT, IsEqualTo(dT_est,'Within',RelativeTolerance(testCase.tol)|AbsoluteTolerance(testCase.tol)), sprintf('Tangent Vector derivative for cubic terrain failed on %s terrain',testCase.labels{n}));
            end
            
        end
        function quinticBasisDerivativeTest(testCase)
            % Import the necessary packages
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.IsEqualTo
            tangentfun = @(w)testCase.outputWrapper2(@(z)testCase.quintic.basis(z),w);
            for n = 1:length(testCase.x)
                % Test point and calculated derivatives
                r = testCase.x{n};
                [~, ~, dN, dT] = testCase.quintic.basis(r);
                
                % Check the derivative of the normal and tangent vectors
                dN_est = testCase.finiteDifference(@(z)testCase.quintic.basis(z), r);
                dT_est = testCase.finiteDifference(@(z)tangentfun(z), r);
                
                
                testCase.verifyThat(dN, IsEqualTo(dN_est,'Within',RelativeTolerance(testCase.tol)|AbsoluteTolerance(testCase.tol)), sprintf('Normal vector derivative for quintic terrain failed on %s terrain', testCase.labels{n}));
                testCase.verifyThat(dT, IsEqualTo(dT_est,'Within',RelativeTolerance(testCase.tol)|AbsoluteTolerance(testCase.tol)), sprintf('Tangent Vector derivative for quintic terrain failed on %s terrain',testCase.labels{n}));
                
            end
        end
    end
    
    
    
    
    methods
        function [g,dg,Hg] = cubic_eval(obj, x)
            [p, D] = obj.cubic.polynomialize(x);
            
            g = obj.cubic.spline_coeff'*p;
            dg = obj.cubic.spline_coeff'*D*p;
            Hg = obj.cubic.spline_coeff'*D^2*p;
            
        end
        function [g, dg, Hg] = quintic_eval(obj,x)
            [p, D] = obj.quintic.polynomialize(x);
            
            g = obj.quintic.spline_coeff'*p;
            dg = obj.quintic.spline_coeff'*D*p;
            Hg = obj.quintic.spline_coeff'*D^2*p;
        end
        function df = finiteDifference(obj, fun, x)
            
            f = fun(x);
            
            df =  cell(1, numel(x));
            % Shifts
            dx = obj.h * eye(numel(x));
            % Central differencing, 6th order accurate
            for n = 1:numel(x)
                f0 = fun(x - 3*dx(:,n));
                f1 = fun(x - 2*dx(:,n));
                f2 = fun(x - dx(:,n));
                f3 = fun(x + dx(:,n));
                f4 = fun(x + 2* dx(:,n));
                f5 = fun(x + 3 * dx(:,n));
                df{n} = (1/60*(f5 - f0) + 3/20* (f1 - f4) + 3/4*(f3 - f2))./obj.h;
            end
            % Collect the outputs together
            s = size(f);
            if s(2) == 1
                dim = 2;
            else
                dim = numel(s) + 1;
            end
            df = cat(dim, df{:});
        end
    end
    methods (Static)
        function f = outputWrapper2(fun, x)
            [~, f] = fun(x);
        end
    end
end

