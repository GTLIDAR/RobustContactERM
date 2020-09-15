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
            if isfield(options, 'contactCostMultiplier') && options.contactCostMultiplier ~= 0
                options.frictionCostMultiplier = options.contactCostMultiplier;
                options.distanceCostMultiplier = options.contactCostMultiplier;
            else
                if ~isfield(options, 'frictionCostMultiplier')
                    options.frictionCostMultiplier = 1;
                end
                if ~isfield(options, 'distanceCostMultiplier')
                    options.distanceCostMultiplier = 1;
                end
            end
            if ~isfield(options, 'uncertainty_source')
                options.uncertainty_source = RobustContactImplicitTrajectoryOptimizer.NO_UNCERTAINTY;
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
            % Check for a distance scaling
            if ~isfield(options, 'distanceScaling')
               options.distanceScaling = 1; 
            end
            % Check for a bias for friction 
            if ~isfield(options, 'ermFrictionBias')
               options.ermFrictionBias = 0; 
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
            
            % Switch out the deterministic solvers for ERM solvers (if
            % necessary)
            switch obj.options.uncertainty_source
                case RobustContactImplicitTrajectoryOptimizer.NO_UNCERTAINTY
                    % Add deterministic costs for friction distance, and
                    % sliding
                    obj = obj.addFrictionConstraint();
                    obj = obj.addSlidingConstraint();
                    obj = obj.addDistanceConstraint();
                  
                case RobustContactImplicitTrajectoryOptimizer.FRICTION_UNCERTAINTY
                    % Add ERM cost for friction, NCC for distance and
                    % sliding
                    obj = obj.addFrictionERM();
                    obj = obj.addSlidingConstraint();
                    obj = obj.addDistanceConstraint();
                    
                case RobustContactImplicitTrajectoryOptimizer.DISTANCE_UNCERTAINTY
                    % Add ERM cost for distance, NCC for friction and
                    % sliding
                    obj = obj.addDistanceERM();
                    obj = obj.addFrictionConstraint();
                    obj = obj.addSlidingConstraint();
     
                case RobustContactImplicitTrajectoryOptimizer.COMBINED_UNCERTAINTY
                    % Add ERM cost for distance and friction, NCC for
                    % sliding
                    obj = obj.addDistanceERM();
                    obj = obj.addFrictionERM();
                    obj = obj.addSlidingConstraint();
                otherwise
                    error('Unknown uncertainty source');
            end
            % Check for relaxed NC Constraints
            if obj.options.nlcc_mode == 5
               % Add in the relaxing variables as a cost
               cost = FunctionHandleObjective(obj.N-1, @(x) obj.relaxed_nc_cost(x));
               cost = cost.setName('RelaxedNLCCost');
               obj = obj.addCost(cost, obj.slack_inds);
            end
        end
    end
    methods
        function cstrVals = calculateContactConstraints(obj, z)
           %% calculateContactConstraints: Helper function for calculating the complementarity constraints
           %
           %    cstrVals = calculateContactConstraints(prob, z) returns a
           %    structure cstrVals containing the constraint defects for
           %    the normal distance, sliding velocity, and friction cone
           %    constraints evaluated at every knot point in the decision
           %    variable list, z
           
           cstrVals = struct();
           cstrVals.normalDistance = zeros(obj.numContacts, obj.N-1);
           cstrVals.slidingVelocity = zeros(obj.numFriction*obj.numContacts, obj.N-1);
           cstrVals.frictionCone = zeros(obj.numContacts, obj.N-1);
           
           for n = 1:obj.N-1
               cstrVals.normalDistance(:,n) = obj.normalDistanceConstraint([z(obj.x_inds(:,n+1)); z(obj.lambda_inds(obj.normal_inds, n))]);
               cstrVals.slidingVelocity(:,n) = obj.slidingVelocityConstraint(z([obj.x_inds(:,n+1); obj.lambda_inds(obj.normal_inds,n); obj.lambda_inds(obj.gamma_inds,n); obj.lambda_inds(obj.tangent_inds,n)]));
               cstrVals.frictionCone(:,n) = obj.frictionConeConstraint(z([obj.lambda_inds(obj.normal_inds,n); obj.lambda_inds(obj.tangent_inds,n); obj.lambda_inds(obj.gamma_inds,n)]));
           end
        end
        %% ---------------- FRICTION CONE -------------------- %%
        function [f, df] = frictionConeDefect(obj, lambda)
            %% FrictionConeDefect: Difference between current and maximum friction force
            %
            %   frictionConeDefect calculates the difference between the
            %   maximum allowed friction force (coeff_friction *
            %   normal_force) and the current value of the friction force.
            %   frictionConeDefect also returns the gradient of the
            %   difference with respect to each of the input force values.
            %
            %   frictionConeDefect is used internally 
            %
            %   Return Values:
            %       f: Nc x 1 vector - friction force difference from maximum, 
            %       df: Nc x Nc matrix -  force difference gradient
            
            % Get the friction coefficient
            nQ = obj.plant.getNumPositions();
            [~, ~, ~, ~, ~, ~, ~, mu] = obj.plant.contactConstraints(zeros(nQ, 1), false, obj.options.active_collision_options);
            
            % Get the forces from lambda
            lambda_N = lambda(obj.normal_inds);
            lambda_T = lambda(obj.tangent_inds);
            % Create a matrix to select and sum up the frictional forces
            % for each contact
            S = zeros(obj.numContacts, obj.numContacts * obj.numFriction);
            for n = 1:obj.numContacts
                S(n,obj.numFriction * (n-1) + 1: obj.numFriction * n) = 1;
            end
            % Calculate the friction cone defect.
            f = mu.*lambda_N - S * lambda_T;
            df = zeros(obj.numContacts, numel(lambda));
            df(:, obj.normal_inds) = diag(mu);
            df(:,obj.tangent_inds) = -S;
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
            
            % Recover the forces from the decision variables
            lambda = zeros(size(y));
            lambda(obj.normal_inds) = y(1:obj.numContacts);
            lambda(obj.tangent_inds)= y(obj.numContacts+1:obj.numContacts*(1+obj.numFriction));
            lambda(obj.gamma_inds) = y(obj.numContacts*(1+obj.numFriction)+1:obj.numContacts*(2 + obj.numFriction));
            % Calculate the defects
            [f, dg] = obj.frictionConeDefect(lambda);
            f = f(:);
            % Re-order the columns of the derivative
            df = zeros(size(dg));
            df(:,1:obj.numContacts) = dg(:,obj.normal_inds);
            df(:,obj.numContacts+1:obj.numContacts*(1 + obj.numFriction)) = dg(:, obj.tangent_inds);
            df(:, obj.numContacts*(1 + obj.numFriction) + 1 : obj.numContacts*(2 + obj.numFriction)) = dg(:, obj.gamma_inds);
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
            gamma = lambda(obj.gamma_inds);
            dgamma = zeros(obj.numContacts, nL);
            dgamma(:,obj.gamma_inds) = eye(obj.numContacts);
            
            % Calculate the ERM Variance
            lambda_N = lambda(obj.normal_inds);
            sigma = obj.options.frictionVariance * lambda_N + obj.options.ermFrictionBias;
            dsigma = zeros(obj.numContacts, nL);
            dsigma(:, obj.normal_inds) = eye(obj.numContacts) * obj.options.frictionVariance;
            
            % Get the ERM cost
            [f, df] = obj.ermCost(gamma, z, sigma);
            % Calculate the total differential cost
            df = df * [dgamma; dz; dsigma];
            % Sum over all contacts
            f = obj.options.frictionCostMultiplier * sum(f, 1);
            df = obj.options.frictionCostMultiplier * sum(df, 1);
            
        end
        %% --------------- NORMAL DISTANCE ----------------- %%
        function [f, df] = normalDistanceDefect(obj, x1)
            %% normalDistanceDefect: normal distance to the terrain at the current state
            %
            %   normalDistanceDefect returns the normal distance to the
            %   terrain, as part of a nonlinear complementarity constraint
            %
            %   Syntax:
            %       [f, df] = obj.normalDistanceDefect(x)
            %
            %   Arguments:
            %       obj: RobustContactImplicitTrajectoryOptimizer instance
            %       x: the state vector of the system
            %
            %   Return value:
            %       f: the normal distance between the terrain and the
            %       contact points
            %       df: the gradient of the normal distance with respect to
            %       the state
            nQ = obj.plant.getNumPositions();
            q = x1(1:nQ);
            [f, ~, ~, ~, ~, ~, ~, ~, df] = obj.plant.contactConstraints(q, false, obj.options.active_collision_options);
            df = [df, zeros(obj.numContacts, nQ)];
        end
        function [f, df] = normalDistanceConstraint(obj, y)
            %% normalDistanceConstraint: complementarity constraint for normal distance and force
            %
            %   normalDistanceConstraint models the complementarity
            %   relationship between the normal force and the shortest
            %   normal distance.
            %
            %   Syntax:
            %       [f, df] = obj.normalDistanceConstraint(y)
            %
            %   Arguments:
            %       obj: RobustContactImplicitTrajectoryOptimizer instance
            %       y: a decision variable list, ordered as [x, lambdaN],
            %       where x is the state and lambdaN is the normal force
            %
            %   Return value:
            %       f: the normal distance between the terrain and the
            %       contact points
            %       df: the gradient of the normal distance with respect to
            %       the decision variables.            
            nX = obj.plant.getNumStates();
            x = y(1:nX);
            [f, df] = obj.normalDistanceDefect(x);
            f = obj.options.distanceScaling .* f(:);
            df = obj.options.distanceScaling .* [df, zeros(obj.numContacts)];
        end
        function [f, df] = normalDistanceERMCost(obj, x, lambda)
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
            
            % Get the ERM variables (possibly scaled) and their derivatives
            
            nX = obj.plant.getNumStates();
            nQ = obj.plant.getNumPositions();
            nL = numel(lambda);
            % Get the normal distance and its gradient
            [phi, ~, ~, ~, ~, ~, ~, ~, Jn] = obj.plant.contactConstraints(x(1:nQ), false, obj.options.active_collision_options);
            % Expand the derivative of the distance function
            dphi = [zeros(obj.numContacts, 1), Jn, zeros(obj.numContacts, nQ + nL)];
            dlambda = [zeros(obj.numContacts,1),zeros(obj.numContacts, nX), eye(obj.numContacts)];
            
            % Calculate the variance of the ERM
            sigma = ones(obj.numContacts, 1) * obj.options.heightVariance;
            dsigma = zeros(obj.numContacts,1 + numel(x) + numel(lambda));
            
            % Scale the distance function
            phi = obj.options.distanceScaling * phi;
            dphi = obj.options.distanceScaling * dphi;
            
            % Calculate the erm Cost
            [f, df] = obj.ermCost(lambda, phi, sigma);
            
            % Calculate the total differential cost
            df = df * [dlambda; dphi; dsigma];           
            % Sum over all contacts
            f = obj.options.distanceCostMultiplier * sum(f, 1);
            df = obj.options.distanceCostMultiplier * sum(df, 1);
            
        end
        %% ---------------- SLIDING VELOCITY --------------- %%
        function [f, df] = slidingVelocityDefect(obj, x1, lambda)
            %% slidingVelocityDefect: Nonlinear Complementarity function modeling sliding velocity constraints
            %
            %  [f, df] = slidingVelocityDefect(obj, x, lambda) calculates
            %  the difference between the velocity slack variable and then
            %  tangential velocity in local contact coordinates.
            %
            %   Arguments:
            %       Obj: A RobustContactImplicitTrajectoryOptimizer
            %       instance
            %       x1: the state variables
            %       lambda: force variables, ordered as [lambdaN_1, lambdaT_1,
            %       gamma_1, lambdaN_2, ...], where lambdaN_i is normal
            %       force at the ith contact point, lambdaT_i is the
            %       friction forces, and gamma_i is the sliding velocity
            %       slack variable.
            %
            %   Return Values:
            %       f: the sliding velocity constraint value
            %       df: the gradient of the sliding velocity constraint
            %       with respect to the state and force variables.
            
            % Get the Tangential Friction Basis
            nX = obj.plant.getNumStates();
            nQ = obj.plant.getNumPositions();
            q = x1(1:nQ);
            dq = x1(nQ+1:end);
            [~, ~, ~, ~, ~, ~, ~, ~,~, Jt, ~, dJt] = obj.plant.contactConstraints(q, false, obj.options.active_collision_options);
            % Reshape the tangential basis so it's easier to use
            % Stack the vectors from different directions together
            Jt = cat(1,Jt{:});
            % Organize the derivatives to match the Jacobian
            for n = 1:length(dJt)
               dJt{n} = reshape(dJt{n}', [nQ, nQ, obj.numContacts]);
               dJt{n} = permute(dJt{n}, [3 1 2]);
            end
            dJt = cat(1,dJt{:});
            % Separate out the forces
            nL = numel(lambda);
            gamma = lambda(obj.gamma_inds);
            % Loop over the number of contacts and calculate the sliding
            % velocity defect
            f = zeros(obj.numContacts * obj.numFriction, 1);
            df = zeros(obj.numContacts * obj.numFriction, nX + nL);
            for n = 1:obj.numContacts
                % Indices for the current force variable
                rangeIdx = obj.numFriction * (n-1) + 1 : (obj.numFriction * n);
                % NCP Variables
                f(rangeIdx) = gamma(n) + Jt(n:obj.numContacts:end,:) * dq;
                dJt_dq = squeeze(sum(dJt(n:obj.numContacts:end,:,:) .* dq', 2));
                df(rangeIdx, 1:nX) = [dJt_dq, Jt(n:obj.numContacts:end,:)];
                df(rangeIdx, nX + obj.gamma_inds(n)) = 1;
            end
        end
        function [f, df] = slidingVelocityConstraint(obj,y)
            %% SlidingVelocityConstraint: Complementarity Constraint on Sliding velocity and friction force
            %
            %   slidingVelocityConstraint implements the complementarity
            %   constraint between the sliding velocity and the friction
            %   force.
            %
            %   Arguments:
            %       obj: A RobustContactImplicitTrajectoryOptimizer
            %       instance
            %       y: A list of decision variables, ordered as [x,
            %       lambdaN, gamma, lambdaT], where x is the state, lambdaN
            %       is the normal force, gamma is the sliding velocity
            %       slack variable, and lambdaT is the friction force.
            %
            %   Return values:
            %       f: The value of the nonlinear complementarity function,
            %       the difference between the sliding velocity slack
            %       variable and the end effector velocity projected onto
            %       the tangent space of the contact surface
            %       df: The gradient of the complementarity function with
            %       respect to the decision variables.
            
            nX = obj.plant.getNumStates();
            x = y(1:nX);
            % Re-order the lambda-variables so we can use them
            lambda = zeros(1,obj.numContacts*(2 + obj.numFriction));
            lambda(obj.normal_inds) = y(nX+1:nX+obj.numContacts);
            lambda(obj.tangent_inds) = y(nX + 2*obj.numContacts + 1:nX + obj.numContacts*(2 + obj.numFriction));
            lambda(obj.gamma_inds) = y(nX + obj.numContacts +1: nX + 2*obj.numContacts);            
            % Get the defects
            [f, df] = obj.slidingVelocityDefect(x, lambda);
            f = f(:);
            % Re-order the columns of df_lambda
            df_lambda = df(:,nX+1:end);
            df_y = zeros(size(df_lambda));
            df_y(:, 1:obj.numContacts) = df_lambda(:, obj.normal_inds);
            df_y(:, obj.numContacts+1 : 2*obj.numContacts) = df_lambda(:, obj.gamma_inds);
            df_y(:, 2*obj.numContacts+1 : obj.numContacts*(2+obj.numFriction)) = df_lambda(:, obj.tangent_inds);
            df(:,nX+1:end) = df_y;
        end
        function [f, df] = slidingVelocityCost(obj, x1, lambda)
            %% SLIDINGVELOCITYMINCOST: Deterministic cost for the sliding velocity constraints
            %
            %   SlidingVeocityCost models the sliding velocity constraints
            %   as a residual cost using the MIN residual function. The
            %   residual is 0 when the constraint is satisfied
            %
            %   Arguments:
            %       OBJ: RobustContactImplicitTrajectoryOptimizer instance
            %       x1:  The state vector at the current knot point
            %       lambda: The force vector at the next knot point,
            %       including the velocity slack variables.
            %
            %   Return Values:
            %       f:  The complementary resiudal function value
            %       df: The gradient of the residual function with
            %       respect to the states and the forces.
            
            [z, dz] = obj.slidingVelocityDefect(x1, lambda);
            nX = obj.plant.getNumStates();
            % Get the tangential forces
            lambda_T = lambda(obj.tangent_inds);
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
    
    methods (Access = protected)
        function obj = addDistanceERM(obj)
            %% addDistanceERM: add the ERM cost for uncertain terrain height to the problem
                        
            grad_level = 1; % Gradients are provided
            % Add ERM cost for distance
            % Create the distance objective
            distance = FunctionHandleObjective(obj.plant.getNumStates() + obj.numContacts, @obj.normalDistanceERMCost, grad_level);
            distance = distance.setName(sprintf('DistanceERMCost'));
            % Add the objective at every point
            for i = 1:obj.N-1
                distanceIdx = {obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds,i)};
                obj = obj.addCost(distance, distanceIdx);
            end
        end
        function obj = addFrictionERM(obj)
            %% addFrictionERM: add the ERM cost for uncertain friction coefficient to the problem
            
            % Add ERM cost for friction, NCP cost for distance
            grad_level = 1;
            nL = obj.numContacts * (2 + obj.numFriction);
            friction = FunctionHandleObjective(nL ,@obj.frictionConeERMCost, grad_level);
            friction = friction.setName(sprintf('FrictionERMCost'));
            % Add the objective at every point
            for i = 1:obj.N-1
                frictionIdx = obj.lambda_inds(:,i);
                obj = obj.addCost(friction, frictionIdx);
            end
            % If desired, also add in the friction constraint
            if obj.options.ermMode == RobustContactImplicitTrajectoryOptimizer.ERM_COMBINED
               obj = obj.addFrictionConstraint(); 
            end
        end
        function obj = addSlidingConstraint(obj)
            %% addSlidingConstraint: adds the complementarity constraint for sliding velocity to the problem
            % Create a Comlementarity Constraint for sliding
            sliding = RelaxedNonlinearComplementarityConstraint(@obj.slidingVelocityConstraint, obj.plant.getNumStates()+2*obj.numContacts, obj.numContacts*obj.numFriction, obj.options.nlcc_mode, obj.options.compl_slack);
            sliding.constraints{1} = sliding.constraints{1}.setName('FrictionForceNonneg');
            sliding.constraints{2} = sliding.constraints{2}.setName('TangentVelocityNonneg');
            sliding.constraints{3} = sliding.constraints{3}.setName('TangentVelocityCompl');
            if sliding.n_slack > 0
               slack_idx = zeros(sliding.n_slack, obj.N-1); 
            end
            % Add the constraint at every knot point
            for i = 1:obj.N-1
                % Order the indices such that the argument is 
                %   [x, lambdaN, gamma, lambdaT]
                if obj.options.nlcc_mode == 5
                    slidingIdx = [obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds, i); obj.lambda_inds(obj.gamma_inds, i); obj.lambda_inds(obj.tangent_inds, i); obj.slack_inds(i)];
                else
                    slidingIdx = [obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds, i); obj.lambda_inds(obj.gamma_inds, i); obj.lambda_inds(obj.tangent_inds, i)];
                end
                if sliding.n_slack > 0
                    slack_idx(:,i) = (obj.num_vars + 1: obj.num_vars + sliding.n_slack)';
                end
                obj = obj.addConstraint(sliding, slidingIdx);
            end
            if sliding.n_slack > 0
                obj.slack_inds = [obj.slack_inds; slack_idx];
            end
        end
        function obj = addFrictionConstraint(obj)
            %% addFrictionConstraint: adds the complementarity constraint for the friction cone to the optimization problem
            % Create a complementarity constraint for the friction cone
             friction = RelaxedNonlinearComplementarityConstraint(@obj.frictionConeConstraint, obj.numContacts*(1+obj.numFriction), obj.numContacts, obj.options.nlcc_mode, obj.options.compl_slack);
             friction.constraints{1} = friction.constraints{1}.setName('SlidingNonneg');
             friction.constraints{2} = friction.constraints{2}.setName('FrictionConeNonneg');
             friction.constraints{3} = friction.constraints{3}.setName('FrictionConeCompl');
             % Add the constraint to every knot point
             if friction.n_slack > 0
                slack_idx = zeros(friction.n_slack, obj.N-1); 
             end
             for i = 1:obj.N-1
                 % Reshape all the lambda_inds so the forces are grouped together as [normal, tangential, slack]
                 if obj.options.nlcc_mode == 5
                     frictionIdx =  [obj.lambda_inds(obj.normal_inds,i); obj.lambda_inds(obj.tangent_inds,i); obj.lambda_inds(obj.gamma_inds,i); obj.slack_inds(i)];
                 else
                     frictionIdx =  [obj.lambda_inds(obj.normal_inds,i); obj.lambda_inds(obj.tangent_inds,i); obj.lambda_inds(obj.gamma_inds,i)];
                 end
                 if friction.n_slack > 0
                    slack_idx(:,i) = (obj.num_vars + 1 : obj.num_vars + friction.n_slack)'; 
                 end
                 obj = obj.addConstraint(friction, frictionIdx);
             end
             if friction.n_slack > 0
                 obj.slack_inds = [obj.slack_inds; slack_idx];
             end
        end
        function obj = addDistanceConstraint(obj)
            %% addDistanceConstraint: adds the complementarity constraint for normal distance to the optimization problem 
            % Create a complementarity constraint for the normal distance
            distance = RelaxedNonlinearComplementarityConstraint(@obj.normalDistanceConstraint, obj.plant.getNumStates(), obj.numContacts, obj.options.nlcc_mode, obj.options.compl_slack);
            distance.constraints{1} = distance.constraints{1}.setName('NormalForceNonNeg');
            distance.constraints{2} = distance.constraints{2}.setName('NormalDistanceNonNeg');
            distance.constraints{3} = distance.constraints{3}.setName('NormalDistanceCompl');
            if distance.n_slack > 0
               slack_idx = zeros(distance.n_slack, obj.N-1); 
            end
            % Add the constraint at every knot point
            for i = 1:obj.N-1
                % For distance, we only need [x, lambdaN];
                if obj.options.nlcc_mode == 5
                    distanceIdx = [obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds,i); obj.slack_inds(i)];
                else
                    distanceIdx = [obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds,i)];
                end
                if distance.n_slack > 0
                   slack_idx(:,i) = (obj.num_vars + 1:obj.num_vars + distance.n_slack)'; 
                end
                obj = obj.addConstraint(distance, distanceIdx);
            end
            if distance.n_slack > 0
                obj.slack_inds = [obj.slack_inds; slack_idx];
            end
        end
    end
    methods (Static)
        %% --------------- Residual Functions -------------------- %%
        function [f, df] = ncpResiduals(x, y, dx, dy)
            %% NCPRESIDUALS: Nonlinear complementarity problem function evaluation
            %
            %   ncpResiduals returns the values of a NCP function and their
            %   gradients. ncpResiduals implements the MIN NCP residual
            %   function.
            %       
            %   Given two variables of the same length x, y, ncpResiduals
            %   returns the elementwise minimum of x and y in f. Given the
            %   gradients of x and y with respect to some third unspecified
            %   parameter, ncpResiduals also returns the gradient of f with
            %   repect to that parameter in df. 
            
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
        function [f, df, Hf] = ermCostGaussian(x, mu, sigma)
            %% ERMCOSTGAUSSIAN: Implementation of the ERM cost function for Gaussian uncertainty
            %
            %   [f, df, Hf] = ermCostGaussian(x,mu, sigma) returns the
            %   Expected Residual for Gaussian Distributed complementarity
            %   variables. In the input, x is the determinsitic variable,
            %   mu is the mean of the Gaussian distribution, and sigma is
            %   the standard deviation of the Gaussian distribution.
            %
            %   ermCostGaussian returns the expected residuals under the
            %   MIN residual function in f, as well as the gradient of f
            %   with respect to all inputs in df, and the associated
            %   Hessian Hf.          
            
            % Column vectorize
            x = x(:);
            mu = mu(:);
            sigma = sigma(:);
            
            nX = numel(x);
            % Initialize the Hessian
            Hf = zeros(nX, 3*nX, 3*nX);
            % Check for degenerate distributions
            degenerate = (sigma <= 0);
            sigma(degenerate) = 0;
            % Initialize the pdf and cdf values
            pdf = zeros(length(x), 1);
            cdf = zeros(length(x), 1);
            % Calculate pdf and cdf for nondegenerate variables
            pdf(~degenerate) = normpdf(x(~degenerate), mu(~degenerate), sigma(~degenerate));
            cdf(~degenerate) = normcdf(x(~degenerate), mu(~degenerate), sigma(~degenerate));
            % Include the limiting case for CDF when x > mu
            cdf(and(degenerate, x > mu)) = 1;
            
            % Calculate the function values
            f = x.^2 - sigma.^2 .* (x + mu) .* pdf + (sigma.^2 + mu.^2 - x.^2) .* cdf;
            % Now the gradients (organized as [df/dx, df/dmu, df/dsigma]);
            df = [diag(2.*x.*(1 - cdf)), diag(2.*(mu.*cdf - sigma.^2 .* pdf)), diag(2.*sigma.*(cdf - x.*pdf))];
            % And the Hessians
            tau = (x - mu)./sigma;
            tau(degenerate) = 0;
            sigma(degenerate) = 1;      % Avoid divide by zeros
            f_xx = 2.*(1 - x.*pdf - cdf);
            f_mm = 2.*(cdf - x.*pdf);
            f_ss = -2 .* ((x-mu) + x.*tau.^2).*pdf + 2 .* cdf;
            f_xm = 2.*x.*pdf;
            f_xs = 2.*x.*tau.*pdf;
            f_ms = -2.*sigma.*(1 + x.*tau./sigma).*pdf; 
            % Deal the values of the hessian to the matrix
            I = eye(nX) .* reshape(eye(nX), nX, 1, nX);
            % Deal in the diagonal values
            Hf(:,1:nX,1:nX) = 0.5 * I .* f_xx;
            Hf(:,nX+1:2*nX, nX+1:2*nX) = 0.5 * I .* f_mm;
            Hf(:,2*nX+1:end, 2*nX+1:end) = 0.5 * I .* f_ss;
            % Deal in the upper triangle
            Hf(:,1:nX, nX+1:2*nX) = I.*f_xm;
            Hf(:,1:nX, 2*nX+1:end) = I.*f_xs;
            Hf(:,nX+1:2*nX, 2*nX+1:end) = I.*f_ms;
            % Fill in the lower triangle
            Hf = Hf + permute(Hf, [1,3,2]);
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