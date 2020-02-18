classdef RobustContactImplicitTrajectoryOptimizer < ContactImplicitTrajectoryOptimizer
    
    properties (Constant)
        % Methods for handling parameteric uncertainty
        NO_UNCERTAINTY = 0;
        FRICTION_UNCERTAINTY = 1;
        DISTANCE_UNCERTAINTY = 2;
        COMBINED_UNCERTAINTY = 3;
        % DISTRIBUTION
        GAUSSIAN = 1;    % Use the Gaussian distribution over the unknowns
        LOGISTIC = 2;    % Use the Logistic distribution over the unknowns
        % RESIDUAL FUNCTION
        NONE = 1;        % Use the traditional nonlinear complementarity solveriplier
        MIN = 2;         % Use a min function as the solver - DO NOT USE
    end
    properties
        ermCost;    % Function handle to switch between distributions
    end
    
    methods
        function obj = RobustContactImplicitTrajectoryOptimizer(plant, N, duration, options)
            %% RobustContactImplicitTrajectoryOptimizer: Constructor method
            %
            %   RobustContactImplicitTrajectoryOptimizer constructs an
            %   instance of the RobustContactImplicitTrajectoryOptimizer
            %   class. The constructor sets some options specific to the
            %   class before passing the bulk of object construction to the
            %   parent class, ContactImplicitTrajectoryOptimization
            
            % Set the options specific to the class
            if nargin < 4
                options = struct();
            end
            if ~isfield(options, 'contactCostMultiplier')
                options.contactCostMultiplier = 1;
            end
            if ~isfield(options, 'uncertainty_source')
                options.uncertainty_source = RobustContactImplicitTrajectoryOptimizer.NO_UNCERTAINTY;
            end
            if ~isfield(options, 'complementarity_solver')
                options.complementarity_solver = RobustContactImplicitTrajectoryOptimizer.NONE;
            end
            % Check for nominal values of the friction and height variances
            if ~isfield(options, 'frictionVariance')
                options.frictionVariance = 1;
            elseif options.frictionVariance < 0
                error('frictionVariance must be a positive scalar');
            end
            if ~isfield(options, 'heightVariance')
                options.heightVariance = 1;
            elseif options.heightVariance < 0
                error('heightVariance must be a positive scalar');
            end
            % Check for the ERM distribution
            if ~isfield(options, 'distribution')
                options.distribution = RobustContactImplicitTrajectoryOptimizer.GAUSSIAN;
            end
            
            % Pass construction to the parent class
            obj = obj@ContactImplicitTrajectoryOptimizer(plant, N, duration, options);
        end
        function obj = addContactConstraints(obj)
            %% ADDCONTACTCONSTRAINTS: Overloads the method in ContactImplicitTrajectoryOptimizer
            %
            %    addContactConstraints adds the relevant constraints/costs
            %    to the nonlinear program. If there is no uncertainty
            %    source, addContactConstraints executes the parent method,
            %    which adds nonlinear complementarity constraints to the
            %    program. If there is an uncertainty source,
            
            % GRAD_LEVEL: Gradient for use in optimization. Options:
            %    -2: Non-differentiable
            %    -1: Unknown (default)
            %     0: No user gradients
            %     1: First derivatives provided
            grad_level = 1;
            
            % Set the handle to the ERM Cost function
            switch obj.options.distribution
                case RobustContactImplicitTrajectoryOptimizer.GAUSSIAN
                    % Set the cost function handle to the Gaussian
                    % distribution case
                    obj.ermCost = @RobustContactImplicitTrajectoryOptimizer.ermCostGaussian;
                case RobustContactImplicitTrajectoryOptimizer.LOGISTIC
                    % Set the cost function handle to the Logistic
                    % distribution case
                    obj.ermCost = @RobustContactImplicitTrajectoryOptimizer.ermCostLogistic;
                otherwise
                    % Unrecognized cost function
                    error('Unrecognized distribution');
            end
            
            nX = obj.plant.getNumStates();
            nL = obj.numContacts * (2 + obj.numFriction);
            
            % Initialize the index sets
            frictionIdx = cell(1, obj.N-1);
            slidingIdx = frictionIdx;
            distanceIdx = frictionIdx;
            
            switch obj.options.complementarity_solver
                case RobustContactImplicitTrajectoryOptimizer.NONE
                    % Add the sliding constraint
                    sliding = NonlinearComplementarityConstraint_original(@obj.slidingVelocityConstraint, nX+2*obj.numContacts, obj.numContacts*obj.numFriction, obj.options.nlcc_mode);
                    addSliding = @(obj, c) addConstraint(obj, sliding, c);
                    
                    % Normal Distance Constraint
                    distance = NonlinearComplementarityConstraint_original(@obj.normalDistanceConstraint, nX, obj.numContacts, obj.options.nlcc_mode);
                    addDistance = @(obj, c) addConstraint(obj, distance, c);
                    
                    % Friction Cone Contraint
                    friction = NonlinearComplementarityConstraint_original(@obj.frictionConeConstraint, obj.numContacts*(1+obj.numFriction), obj.numContacts, obj.options.nlcc_mode);
                    addFriction = @(obj, c) addConstraint(obj, friction, c);
                    
                    % Create the index sets
                    for i = 1:obj.N-1
                        % Reshape all the lambda_inds so the forces are grouped together as [normal, tangential, slack]
                        lam_ind = obj.force_converter'*obj.lambda_inds(:,i);
                        frictionIdx{i} =  lam_ind;
                        slidingIdx{i} = [obj.x_inds(:,i+1); lam_ind(1:obj.numContacts); lam_ind(obj.numContacts*(1+obj.numFriction)+1:end); lam_ind(obj.numContacts+1:obj.numContacts*(1+obj.numFriction))];
                        distanceIdx{i} = [obj.x_inds(:,i+1); lam_ind(1:obj.numContacts)];
                    end
                    
                case RobustContactImplicitTrajectoryOptimizer.MIN
                    %NOTE: MIN Currently does not return feasible solutions
                    
                    %Sliding constraint
                    sliding = FunctionHandleObjective(nX+nL, @obj.slidingVelocityCost, grad_level);
                    addSliding = @(obj, c) addCost(obj, sliding, c);
                    
                    % Normal Distance Constraint
                    distance = FunctionHandleObjective(nX+obj.numContacts, @obj.normalDistanceCost, grad_level);
                    addDistance = @(obj, c) addCost(obj, distance, c);
                    
                    % Friction Cone Constriant
                    friction = FunctionHandleObjective(nX + nL, @obj.frictionConeCost, grad_level);
                    addFriction = @(obj, c) addCost(obj, friction, c);
                    
                    % Create the index sets
                    for i = 1:obj.N-1
                        frictionIdx{i} = {obj.lambda_inds(:,i)};
                        slidingIdx{i} = {obj.x_inds(:,i+1); obj.lambda_inds(:,i)};
                        distanceIdx{i} = {obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds,i)};
                    end
                    
                otherwise
                    error('Unrecognized complementarity solver option');
            end
            
            % Switch out the deterministic solvers for ERM solvers (if
            % necessary)
            switch obj.options.uncertainty_source
                case RobustContactImplicitTrajectoryOptimizer.NO_UNCERTAINTY
                    % Add deterministic costs for friction and distance - Do
                    % Nothing
                    
                case RobustContactImplicitTrajectoryOptimizer.FRICTION_UNCERTAINTY
                    % Add ERM cost for friction, NCP cost for distance
                    friction = FunctionHandleObjective(nL ,@obj.frictionConeERMCost, grad_level);
                    addFriction = @(obj, c) addCost(obj, friction, c);
                    % Create the index sets
                    for i = 1:obj.N-1
                        frictionIdx{i} = {obj.lambda_inds(:,i)};
                    end
                    
                    
                case RobustContactImplicitTrajectoryOptimizer.DISTANCE_UNCERTAINTY
                    % Add ERM cost for distance, NCP cost for friction
                    distance = FunctionHandleObjective(nX + obj.numContacts, @obj.normalDistanceERMCost, grad_level);
                    addDistance = @(obj, c) addCost(obj, distance, c);
                    % Create the index sets
                    for i = 1:obj.N-1
                        distanceIdx{i} = {obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds,i)};
                    end
                case RobustContactImplicitTrajectoryOptimizer.COMBINED_UNCERTAINTY
                    % Add ERM cost for distance and friction
                    friction = FunctionHandleObjective(nL, @obj.frictionConeERMCost, grad_level);
                    distance = FunctionHandleObjective(nX + obj.numContacts, @obj.normalDistanceERMCost, grad_level);
                    addFriction = @(obj, c) addCost(obj, friction, c);
                    addDistance = @(obj, c) addCost(obj, distance, c);
                    % Create the index sets
                    for i = 1:obj.N-1
                        frictionIdx{i} = {obj.lambda_inds(:,i)};
                        distanceIdx{i} = {obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds,i)};
                    end
                otherwise
                    error('Unknown uncertainty source');
            end
            
            % Add all the costs/constraints to the problem
            for i = 1:obj.N-1
                obj = addFriction(obj, frictionIdx{i});
                obj = addDistance(obj, distanceIdx{i});
                obj = addSliding(obj, slidingIdx{i});
            end
        end
    end
    methods
        %% ---------------- FRICTION CONE -------------------- %%
        function [f, df] = frictionConeDefect(obj, lambda)
            
            % Get the friction coefficient
            nQ = obj.plant.getNumPositions();
            [~, ~, ~, ~, ~, ~, ~, mu] = obj.plant.contactConstraints(zeros(nQ, 1), false, obj.options.active_collision_options);
            
            % Get the forces from lambda
            nL = numel(lambda);
            skip = 2 + obj.numFriction;
            T_idx = 1:nL;
            % Normal force
            N_idx = 1:skip:nL;
            lambda_N = lambda(N_idx);
            
            % Tangential force
            G_idx = skip:skip:nL;
            T_idx([N_idx,G_idx]) = [];
            lambda_T = lambda(T_idx);
            % Create a matrix to select and sum up the frictional forces
            % for each contact
            S = zeros(obj.numContacts, obj.numContacts * obj.numFriction);
            for n = 1:obj.numContacts
                S(n,obj.numFriction * (n-1) + 1: obj.numFriction * n) = 1;
            end
            % Calculate the friction cone defect.
            f = mu*lambda_N - S * lambda_T;
            df = zeros(obj.numContacts, nL);
            df(:, N_idx) = eye(obj.numContacts) * mu;
            df(:,T_idx) = -S;
        end
        function [f, df] = frictionConeConstraint(obj, y)
            % For the friction cone constraint, the variable y is stored as:
            %
            %    y = [lambdaN, lambdaT, gamma], and the decision variable is
            %    gamma (lambdaN and lambdaT are parameters). Thus, the
            %    normal forces for all contacts must be stored together
            %    (lambdaN = [lambdaN_1, lambdaN_2,... lambdaN_M]), and
            %    likewise for the tangential forces and the velocity slack
            %    variables.
            %
            %    This function is intended to be used with the
            %    NonlinearComplementarityConstraint class.
            
            % Re-order the decision variables to be used with the defect
            % function (use a matrix so it's invertible).
            lambda = obj.force_converter*y;
            [f, df] = obj.frictionConeDefect(lambda);
            % Re-order the columns of df
            df = df*obj.force_converter;
        end
        function [f, df] = frictionConeCost(obj, lambda)
            %% FrictionConeCost: Cost Function for enforcing deterministic friction cone constraints
            %
            %
            
            [z, dz] = obj.frictionConeDefect(lambda);
            
            % Sliding velocity slack
            nL = numel(lambda);
            skip = 2 + obj.numFriction;
            G_idx = skip:skip:nL;
            gamma = lambda(G_idx);
            % Get the NCP function variables and their gradients
            dgamma = zeros(obj.numContacts, nL);
            dgamma(:,G_idx) = eye(obj.numContacts);
            % Get the NCP function residuals and their gradient
            [f, df] = obj.ncpResiduals(gamma, z, dgamma, dz);
            % Sum the cost and gradients over all contacts
            f = obj.options.contactCostMultiplier * sum(f, 1);
            df = obj.options.contactCostMultiplier * sum(df, 1);
        end
        function [f, df] = frictionConeERMCost(obj, lambda)
            %% FRICTIONCONEERMCOST: Expected Residual for uncertain friction cones
            %
            %   frictionConeERMCost returns the expected residual for the
            %   ERM problem of a uncertain friction cone, where the
            %   friction coefficient is normally distributed.
            %
            %   [f, df] = frictionConeERMCost(obj, lambda)
            %
            %   where lambda are the force decision variables, stored as:
            %   lambda = [lambda_N1, lambda_T1, gamma_1, ..., gamma_M]
            %
            %       lambda_Nm is the normal force for the mth contact
            %       lambda_Tm are the tangential force components for the
            %       mth contact
            %       gamma_m is the slack variable for sliding velocity for
            %       the mth contact
            
            % Get the friction cone defect and it's gradient
            [z, dz] = obj.frictionConeDefect(lambda);
            % Get the sliding velocity and it's gradient
            nL = numel(lambda);
            skip = 2 + obj.numFriction;
            G_idx = skip:skip:nL;
            gamma = lambda(G_idx);
            dgamma = zeros(obj.numContacts, nL);
            dgamma(:,G_idx) = eye(obj.numContacts);
            
            % Calculate the ERM Variance
            N_idx = 1:skip:nL;
            lambda_N = lambda(N_idx);
            sigma = obj.options.frictionVariance * lambda_N;
            dsigma = zeros(obj.numContacts, nL);
            dsigma(:, N_idx) = eye(obj.numContacts) * obj.options.frictionVariance;
            
            % Get the ERM cost
            [f, df_gamma, df_z, df_sigma] = obj.ermCost(gamma, z, sigma);
            % Calculate the total differential cost
            df = df_gamma * dgamma + df_z * dz + df_sigma * dsigma;
            % Sum over all contacts
            f = obj.options.contactCostMultiplier * sum(f, 1);
            df = obj.options.contactCostMultiplier * sum(df, 1);
            
        end
        %% --------------- NORMAL DISTANCE ----------------- %%
        function [f, df] = normalDistanceDefect(obj, x1)
            % Get the normal distance and its gradient
            nQ = obj.plant.getNumPositions();
            q = x1(1:nQ);
            [f, ~, ~, ~, ~, ~, ~, ~, df] = obj.plant.contactConstraints(q, false, obj.options.active_collision_options);
            df = [df, zeros(obj.numContacts, nQ)];
        end
        function [f, df] = normalDistanceConstraint(obj, y)
            % Here the y variables must be ordered [x, lambdaN]
            %
            %   This function is intended for use with the Nonlinear
            %   ComplementarityConstraint class
            
            nX = obj.plant.getNumStates();
            x = y(1:nX);
            [f, df] = obj.normalDistanceDefect(x);
            df = [df, zeros(obj.numContacts)];
        end
        function [f, df] = normalDistanceCost(obj, x1, lambda)
            %% NormalDistanceCost: Objective function for the nonpenetration constraint
            %
            %   [f, df] = normalDistanceCost(obj, x, lambda) calculates
            %   the objective f for the nonpenetration constraint, given
            %   the current state x1 and the normal force LAMBDA.
            %
            %   normalDistanceCost uses a residual function to calculate
            %   the objective f. The residual function is specified such
            %   that it attains a root when the complementarity conditions
            %   are satisfied:
            %
            %       f(a, b) = 0 iff  a >= 0, b>= 0, ab = 0
            %
            %   the function also returns the derivatives of the objective
            %   with respect to the state X and the normal force LAMBDA:
            %       df = [df/dx, df/dlambda]
            %
            %   Here, we only need the normal force variables, i.e.
            %   lambda = lambda_N
            
            % Get the normal distance and its gradient
            [z, dz] = obj.normalDistanceDefect(x1);
            % Expand to include derivatives wrt lambda
            nQ = obj.plant.getNumPositions();
            nX = obj.plant.getNumStates();
            nL = numel(lambda);
            dz = [dz, zeros(obj.numContacts, nQ + nL)];
            
            % Calculate the derivatives
            dlambda = zeros(obj.numContacts, nX + nL);
            dlambda(:, nX +1:end) = eye(obj.numContacts);
            
            % Use the NCP function to get the elementwise residuals
            [f, df] = obj.ncpResiduals(z, lambda, dz, dlambda);
            
            % Sum the residuals and gradients together for the final cost
            f = obj.options.contactCostMultiplier * sum(f, 1);
            df = obj.options.contactCostMultiplier * sum(df, 1);
        end
        function [f, df] = normalDistanceERMCost(obj, x1, lambda)
            %% NormalDistanceERMCost: Objective function for the nonpenetration constraint under the ERM model
            %
            %   [f, df] = normalDistanceERMCost(obj, x, lambda) calculates
            %   the Expected Residual Minimization (ERM) objective f for
            %   the nonpenetration constraint, given the current state x1
            %   and the normal force LAMBDA.
            %
            %   the function also returns the derivatives of the objective
            %   with respect to the state X and the normal force LAMBDA:
            %       df = [df/dx, df/dlambda]
            %
            %   Here, we only need the normal force variables, i.e.
            %   lambda = lambda_N
                       
            % Get the normal distance and its gradient
            [z, dz] = obj.normalDistanceDefect(x1);
            
            % Expand to include derivatives wrt lambda
            nX = obj.plant.getNumStates();
            nL = numel(lambda);
            dz = [dz, zeros(obj.numContacts, nL)];         
            
            % Get the normal force
            dlambda = zeros(obj.numContacts, nX + nL);
            dlambda(:, nX + 1:end) =  eye(obj.numContacts);
            
            % Calculate the variance of the ERM
            sigma = ones(obj.numContacts, 1) * obj.options.heightVariance;
            dsigma = zeros(obj.numContacts, nX + nL);
            
            % Calculate the erm Cost
            [f, df_lambda, df_z, df_sigma] = obj.ermCost(lambda, z, sigma);
            
            % Calculate the total differential cost
            df = df_lambda * dlambda + df_z * dz + df_sigma * dsigma;
            % Sum over all contacts
            f = obj.options.contactCostMultiplier * sum(f, 1);
            df = obj.options.contactCostMultiplier * sum(df, 1);
            
        end
        %% ---------------- SLIDING VELOCITY --------------- %%
        function [f, df] = slidingVelocityDefect(obj, x1, lambda)
            
            
            % Get the Tangential Friction Basisiplier
            nX = obj.plant.getNumStates();
            nQ = obj.plant.getNumPositions();
            q = x1(1:nQ);
            dq = x1(nQ+1:end);
            [~, ~, ~, ~, ~, ~, ~, ~,~, Jt, ~, dJt] = obj.plant.contactConstraints(q, false, obj.options.active_collision_options);
            % Reshape the tangential basis so it's easier to use
            Jt = cat(3,Jt{:});
            Jt = permute(Jt, [3,2,1]);
            dJt = cat(3,dJt{:});
            dJt = reshape(dJt, [obj.numContacts, nQ, nQ, size(dJt, 3)]);
            dJt = permute(dJt, [4,2,3,1]);
            % Separate out the forces
            nL = numel(lambda);
            skip = 2 + obj.numFriction;
            % Sliding velocity slack
            G_idx = skip:skip:nL;
            gamma = lambda(G_idx);
            % Loop over the number of contacts and calculate the sliding
            % velocity defect
            f = zeros(obj.numContacts * obj.numFriction, 1);
            df = zeros(obj.numContacts * obj.numFriction, nX + nL);
            for n = 1:obj.numContacts
                % Indices for the current force variablesiplier
                rangeIdx = obj.numFriction * (n-1) + 1 : obj.numFriction * n;
                % NCP Variables
                f(rangeIdx) = gamma(n) + Jt(:,:,n) * dq;
                dJt_dq = squeeze(sum(dJt(:,:,:,n) .* dq', 2));
                df(rangeIdx, 1:nX) = [dJt_dq, Jt(:,:,n)];
                df(rangeIdx, nX + G_idx(n)) = 1;
            end
        end
        function [f, df] = slidingVelocityConstraint(obj,y)
            % Here, the y variables should be ordered [x, lambdaN, gamma, lambdaT]
            
            nX = obj.plant.getNumStates();
            x = y(1:nX);
            % Re-order the lambda-variables so we can use them
            lambda = y(nX+1:nX + obj.numContacts*(2+obj.numFriction));
            lambda = [lambda(1:obj.numContacts); lambda(2*obj.numContacts+1:end); lambda(obj.numContacts + 1:2*obj.numContacts)];
            
            lambda = obj.force_converter * lambda;
            
            % Get the defects
            [f, df] = obj.slidingVelocityDefect(x, lambda);
            % Re-order the columns of df_lambda
            df_lambda = df(:,nX+1:end);
            df_lambda = df_lambda * obj.force_converter;
            df_lambda = [df_lambda(:, 1:obj.numContacts), df_lambda(:, obj.numContacts*(obj.numFriction+1)+1:end),df_lambda(:, obj.numContacts+1:obj.numContacts*(obj.numFriction+1))];
            df(:,nX+1:end) = df_lambda;
        end
        function [f, df] = slidingVelocityCost(obj, x1, lambda)
            %% SLIDINGVELOCITYMINCOST: Deterministic cost for the sliding velocity constraints
            
            [z, dz] = obj.slidingVelocityDefect(x1, lambda);
            nX = obj.plant.getNumStates();
            % Get the tangential forces
            nL = numel(lambda);
            skip = 2 + obj.numFriction;
            T_idx = 1:nL;
            N_idx = 1:skip:nL;
            G_idx = skip:skip:nL;
            T_idx([N_idx,G_idx]) = [];
            lambda_T = lambda(T_idx);
            % Gradient of the tangential forces
            dlambda_T = zeros(size(dz));
            dlambda_T(:, nX + T_idx) = eye(length(lambda_T));
            % Get the NCP residuals
            [f, df] = obj.ncpResiduals(lambda_T, z, dlambda_T, dz);
            % Sum the costs and the gradients
            f = obj.options.contactCostMultiplier * sum(f, 1);
            df = obj.options.contactCostMultiplier * sum(df, 1);
        end
    end
    
    methods
        function obj = enableDisplayFunction(obj)
            %% ENABLEDISPLAYFUNCTION: Enables the native display function
            clear('obj.printCostFunction'); %Reset the function and its persistent variables
            obj = obj.addDisplayFunction(@(h, x, u, l)obj.printCostFunction(h, x, u, l));
        end
        function obj = printCostFunction(obj, h, x, u, l)
            %% PRINTCOSTFUNCTIONS: Prints the cost functions at each iteration of the optimization
            %
            %   Note: the display function gets called on every iteration
            %   of the optimization, so we can use persistent variables to
            %   keep track of the iterations
            persistent iteration
            if isempty(iteration)
                iteration = 1;
                fprintf('\n');
            end
            [runningCost, dynamicCstr, distanceCost, frictionCost, slidingCost] = obj.calculateCosts(h, x, t, u, l);       
            % Print the values to the screen
            if rem(iteration, 100) == 1
                % Every 50 iterations, re-print the titlesiplier
                fprintf(' Iteration \t RunningCost \t DynamicCstr \t FrictionCost \t DistanceCost \t SlidingCstr\n');
            end
            fprintf(' %12f \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n', iteration, runningCost, dynamicCstr, frictionCost, distanceCost, slidingCost);
            %Increment the iteration count
            iteration = iteration + 1;
        end
        function [runningCost, dynamicCstr, distanceCost, frictionCost, slidingCost] = calculateCosts(obj, h, x, u, l)
            % Calculate the total ERM Cost
            frictionCost = 0;
            distanceCost = 0;
            slidingCost = 0;
            switch obj.options.uncertainty_source
                case RobustContactImplicitTrajectoryOptimizer.FRICTION_UNCERTAINTY
                    for n = 1:obj.N-1
                        frictionCost = frictionCost + obj.frictionConeERMCost(l(:,n));
                        distanceCost = distanceCost + norm(obj.normalDistanceConstraint([x(:,n+1), l(obj.normal_inds, n)]));
                        slidingCost = slidingCost + norm(obj.slidingVelocityConstraint([x(:,n+1); l(obj.normal_inds,n); l(obj.gamma_inds, n); l(obj.tangent_inds, n)]));
                    end
                case RobustContactImplicitTrajectoryOptimizer.DISTANCE_UNCERTAINTY
                    for n = 1:obj.N-1
                        distanceCost = distanceCost + obj.normalDistanceERMCost(x(:,n+1), l(obj.normal_inds, n));
                        frictionCost = frictionCost + norm(obj.frictionConeConstraint([l(obj.normal_inds, n); l(obj.tangent_inds,n); l(obj.gamma_inds, n)]));
                        slidingCost = slidingCost + norm(obj.slidingVelocityConstraint([x(:,n+1); l(obj.normal_inds,n); l(obj.gamma_inds, n); l(obj.tangent_inds, n)]));
                    end
                case RobustContactImplicitTrajectoryOptimizer.COMBINED_UNCERTAINTY
                    for n = 1:obj.N-1
                        frictionCost = frictionCost + obj.frictionConeERMCost(l(:,n));
                        distanceCost = distanceCost + obj.normalDistanceERMCost(x(:,n+1), l(obj.normal_inds, n));
                        slidingCost = slidingCost + norm(obj.slidingVelocityConstraint([x(:,n+1); l(obj.normal_inds,n); l(obj.gamma_inds, n); l(obj.tangent_inds, n)]));
                    end
                otherwise
                    for n = 1:obj.N-1
                        frictionCost = frictionCost + norm(obj.frictionConeConstraint([l(obj.normal_inds, n); l(obj.tangent_inds,n); l(obj.gamma_inds, n)]));
                        distanceCost = distanceCost + norm(obj.normalDistanceConstraint([x(:,n+1), l(obj.normal_inds, n)]));
                        slidingCost = slidingCost + norm(obj.slidingVelocityConstraint([x(:,n+1); l(obj.normal_inds,n); l(obj.gamma_inds, n); l(obj.tangent_inds, n)]));
                    end
            end
            % Calculate the total cost functional (control objective)
            runningCost = 0;
            dynamicCstr = 0;
            switch obj.options.integration_method
                case RobustContactImplicitTrajectoryOptimizer.FORWARD_EULER
                    for n = 1:obj.N-1
                       runningCost = runningCost + h(n)*obj.cost_handle(h(n),x(:,n),u(:,n)); 
                       dynamicCstr = dynamicCstr + norm(obj.forward_constraint_fun(h(n), x(:,n), x(:,n+1), u(:,n), l(obj.force_inds,n)));
                    end
                case RobustContactImplicitTrajectoryOptimizer.BACKWARD_EULER
                    for n = 1:obj.N-1
                       runningCost = runningCost + h(n)*obj.cost_handle(h(n), x(:, n+1), u(:, n)); 
                       dynamicCstr = dynamicCstr + norm(obj.forward_constraint_fun(h(n), x(:,n), x(:,n+1), u(:,n), l(obj.force_inds,n)));
                    end
                case RobustContactImplicitTrajectoryOptimizer.MIDPOINT
                    for n = 1:obj.N-1
                        runningCost = runningCost + h(n)*obj.midpoint_running_fun(@obj.cost_handle, h(n), x(:,n), x(:,n+1), u(:,n), u(:, n+1));
                        dynamicCstr = dynamicCstr + norm(obj.forward_constraint_fun(h(n), x(:,n), x(:,n+1), u(:,n), u(:,n+1), l(obj.force_inds,n)));
                    end
                case RobustContactImplicitTrajectoryOptimizer.SEMI_IMPLICIT
                    for n = 1:obj.N-1
                        runningCost = runningCost + h(n)*obj.cost_handle(h(n),x(:,n),u(:,n));
                        dynamicCstr = dynamicCstr + norm(obj.forward_constraint_fun(h(n), x(:,n), x(:,n+1), u(:,n), l(obj.force_inds,n)));
                    end
            end     
        end
    end
    
    
    methods (Static)
        function [f, df] = ncpResiduals(x, y, dx, dy)
            %% NCPRESIDUALS: Nonlinear complementarity problem function evaluation
            %
            %   ncpResiduals returns the values of a NCP function and their
            %   gradients.
            
            % Determine which variables are smaller
            xmin_idx = x < y;
            % Cost function: min-squared
            f = zeros(numel(x),1);
            df = zeros(numel(x), size(dx, 2));
            if any(xmin_idx)
                f(xmin_idx) = x(xmin_idx).^2;
                df(xmin_idx,:) =  2.*x(xmin_idx) .* dx(xmin_idx,:);
            end
            if any(~xmin_idx)
                f(~xmin_idx) = y(~xmin_idx).^2;
                df(~xmin_idx,:) = 2.*y(~xmin_idx) .*dy (~xmin_idx,:);
            end
        end
        function [f, dfx, dfmu, dfsigma] = ermCostGaussian(x, mu, sigma)
            %% ERMCOSTGassian: Cost function for the Expected Residual Minimization with Gaussian variables
            %
            %
            
            % Initialize the outputs
            f = zeros(length(x),1);
            dfx = zeros(length(x),1);
            dfmu = dfx;
            dfsigma = dfx;
            % Filter out any degenerate distributions (no variance cases)
            degenerate = (sigma == 0);
            % Save the degenerate means for the limiting case
            mu_degenerate = mu(degenerate);
            x_degenerate = x(degenerate);
            % Filter out the remaining variables
            
            x = x(~degenerate);
            sigma = sigma(~degenerate);
            mu = mu(~degenerate);
            if any(~degenerate)
                % Calculate the pdf and cdf values
                pdf = normpdf(x, mu, sigma);
                cdf = normcdf(x, mu, sigma);
                % The ERM cost function for Gaussian variables with nonzero
                % variance
                f(~degenerate) = x.^2 - sigma.^2 .* (x + mu) .* pdf + (sigma.^2 + mu.^2 - x.^2) .* cdf;
                % Calculate the derivaties of pdf and cdf wrt x, mu, and sigma
                tau = (x - mu)./sigma;
                % The derivatives of pdf
                dp_sigma = 1./sigma .* pdf .*(tau.^2 -1);
                dp_mu = tau.*pdf ./ sigma;
                dp_x = - dp_mu;
                % The derivatives of cdf
                dc_sigma = -pdf .* tau;
                dc_mu = -pdf;
                dc_x = pdf;
                % The derivatives
                dfx(~degenerate) = 2*x - sigma.^2 .* (pdf + (x + mu) .* dp_x) - 2.*x .* cdf + (sigma.^2 + mu.^2 - x.^2) .* dc_x;
                dfmu(~degenerate) = -sigma.^2 .* (pdf + (x + mu).*dp_mu) + 2.*mu .* cdf + (sigma.^2 + mu.^2 - x.^2 ).* dc_mu;
                dfsigma(~degenerate) = -2.*sigma .* (x + mu) .* pdf - sigma.^2 .* (x + mu) .* dp_sigma + 2.*sigma .* cdf + (sigma.^2 + mu.^2 - x.^2) .* dc_sigma;
            end
            if any (degenerate)
                % Handle the degenerate case
                x_idx = (x_degenerate < mu_degenerate);
                g = zeros(length(x_degenerate), 1);
                
                dg_x = g;
                dg_mu = g;
                
                g(x_idx) = x_degenerate.^2;
                g(~x_idx) = mu_degenerate.^2;
                
                dg_x(x_idx) = 2.*x_degenerate;
                dg_mu(x_idx) = 2.*mu_degenerate;
                
                f(degenerate) = g;
                dfx(degenerate) = dg_x;
                dfmu(degenerate) = dg_mu;
            end
            dfx = diag(dfx);
            dfmu = diag(dfmu);
            dfsigma = diag(dfsigma);
        end     
        function [f, dfx, dfmu, dfsigma] = ermCostLogistic(x, mu, sigma)
            %% ERMCostLogistic: Expected Residual Cost for Logistic Random Variables
            
            % Check for the degenerate case (sigma = 0)
            degenerate = (sigma == 0);
            x_degenerate = x(degenerate);
            mu_degenerate = mu(degenerate);
            % Set up the return values
            f = zeros(length(x), 1);
            dfx = f;
            dfmu = f;
            dfsigma = f;
            % Filter out the degenerate variables
            x = x(~degenerate);
            mu = mu(~degenerate);
            sigma = sigma(~degenerate);
            if any(~degenerate)
                % Calculate the ERM Cost for Logistically Distributed variables
                z = exp((x - mu)./sigma);
                Li2 = -dilog(1 + z);        % The Dilogarithm function
                
                % The cost function
                f(~degenerate) = x.^2 - 2.*x.*sigma.*log(1+z) + 2.*sigma.^2 .* Li2;
                
                % The gradients
                dfx(~degenerate) = 2.*x./(1+z);
                dfmu(~degenerate) = 2.*x .* z ./(1+z) - 2.*sigma.* log(1+z);
                dfsigma(~degenerate) = (2*mu - 4 *x) .* log(1+z) + 2*x .* (x-mu)./ sigma .* (z ./(1+z)) + 4 * sigma .* Li2;
            end
            % In the limiting case, the logistic ERM reverts back to the
            % min function
            if any(degenerate)
                g = zeros(length(x_degenerate),1);
                dg_x = g;
                dg_mu = g;
                xidx = x_degenerate < mu_degenerate;
                
                g(xidx) = x_degenerate.^2;
                g(~xidx) = mu_degenerate.^2;
                
                dg_x(xidx) = 2.*x_degenerate;
                dg_mu(~xidx) = 2.*mu_degenerate;
                
                f(degenerate) = g;
                dfx(degenerate) = dg_x;
                dfmu(degenerate) = dg_mu;
            end
            % Convert the gradients into matrices
            dfx = diag(dfx);
            dfmu = diag(dfmu);
            dfsigma = diag(dfsigma);
        end
    end
end