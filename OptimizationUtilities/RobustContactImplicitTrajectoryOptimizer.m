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
        % SCALING
        NOSCALE = 0;
        CONSTANT = 1;
        INERTIA = 2;
        % ERM MODE
        ERM_OBJECTIVE = 1;
        ERM_SEPARATED = 2;
        ERM_CONSTRAINT = 3;
        ERM_COMBINED = 4;
        % ERM SPACE SCALING
        ERM_LINEARSPACE = 1;
        ERM_LOGSPACE = 2;
    end
    properties
        ermCost;    % Function handle to switch between distributions
        distance_slacks;
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
            % Check for ERM Scaling
            if ~isfield(options, 'ermScaling')
               options.ermScaling = RobustContactImplicitTrajectoryOptimizer.NOSCALE; 
            end
            % Check for the ERM Cost Method
            if ~isfield(options, 'ermMode')
                options.ermMode = RobustContactImplicitTrajectoryOptimizer.ERM_COMBINED;
            end
            % Check for a distance scaling
            if ~isfield(options, 'distanceScaling')
               options.distanceScaling = 1; 
            end
            % Check for an ERM space scaling
            if ~isfield(options, 'ermSpace')
               options.ermSpace = RobustContactImplicitTrajectoryOptimizer.ERM_LINEARSPACE; 
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
               obj = obj.addCost(FunctionHandleObjective(obj.N-1, @(x) obj.relaxed_nc_cost(x)), obj.relax_inds);
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
            f = f(:);
            % Re-order the columns of df
            df = df*obj.force_converter;
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
            [f, df] = obj.ermCost(gamma, z, sigma);
            % Calculate the total differential cost
            df = df * [dgamma; dz; dsigma];
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
            f = obj.options.distanceScaling .* f(:);
            df = obj.options.distanceScaling .* [df, zeros(obj.numContacts)];
        end
        function [phi, lambda_s, dphi, dlambda_s, d2phi, d2lambda_s] = getNormalDistanceERMVariables(obj, h, x, lambda)
            %% getNormalDistanceERMVariables: shared function for calculating the ERM Variables (with scaling)
            
            nX = obj.plant.getNumStates();
            nQ = obj.plant.getNumPositions();
            nL = numel(lambda);
            % Get the normal distance and its gradient
            [phi, ~, ~, ~, ~, ~, ~, ~, Jn, ~, dJn] = obj.plant.contactConstraints(x(1:nQ), false, obj.options.active_collision_options);
            % Expand the derivative of the distance function
            dphi = [zeros(obj.numContacts, 1), Jn, zeros(obj.numContacts, nQ + nL)];
            % Expand the hessian of the distance
            d2phi = zeros(obj.numContacts, 1 + 2*nQ + nL, 1 + 2 * nQ + nL);
            d2phi(:,2:nQ+1, 2:nQ+1) = dJn;
            % Initialize the Hessian for the normal force
            d2lambda_s = zeros(obj.numContacts, 1 + nX + nL, 1 + nX + nL);
            % Apply Scaling, if desired
            switch obj.options.ermScaling
                case RobustContactImplicitTrajectoryOptimizer.CONSTANT
                    % Scale the contact force by the total mass and the timestep
                    m = obj.plant.totalMass();
                    lambda_s = h^2 * lambda/m;
                    dlambda_s = [2*h * lambda / m, zeros(obj.numContacts, nX), h^2/m*ones(nL,1)];
                    d2lambda_s(:,:,1) = [2/m * lambda, zeros(obj.numContacts,nX), 2*h/m * ones(nL, 1)];
                    d2lambda_s(:,1,2+nX:end) = 2/m * eye(nL);
                case RobustContactImplicitTrajectoryOptimizer.INERTIA
                    % Calculate the effective inertia
                    [M, ~, ~, dM] = obj.plant.manipulatorDynamics(x(1:nQ), x(nQ+1:end));
                    R = chol(M);
                    Jr = Jn/R;
                    % Scale the force
                    S = Jr * Jr';                   % Mass scale matrix (inverse matrix in normal coordinates)
                    lambda_s = h.^2 * S * lambda;   % Eliminate mass and time units
                    
                    % Calculate the derivatives of the scaling matrix
                    dM = reshape(dM', nQ*[1,1,1]);
                    dJn = reshape(dJn', [obj.numContacts, nQ, nQ]);
                    diM = zeros(nQ*[1, 1, 1]);
                    dS = zeros(nL, nL, nQ);
                    for n = 1:nQ
                        diM(:,:,n) = -((R\(R'\dM(:,:,n)))/R)/R';
                        dS(:,:,n) = dJn(:,:,n) * (R\Jr') + Jn * diM(:,:,n) * Jn' + (Jr/R') * dJn(:,:,n)';
                    end
                    dS_l = sum(dS .* lambda', 2);
                    dS_l = reshape(dS_l, size(dS_l, 1), size(dS_l, 3));
                    
                    dlambda_s = [2 * h * S * lambda, h^2 * dS_l, zeros(obj.numContacts, nQ), h^2 * S];
                    
                otherwise
                    % No scaling
                    lambda_s = lambda;
                    dlambda_s = [zeros(obj.numContacts,1),zeros(obj.numContacts, nX), eye(obj.numContacts)];
                    
            end
            
        end
        function [f, df] = normalDistanceERMCost(obj, h, x, lambda)
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
            [phi, lambda_s, dphi, dlambda_s] = obj.getNormalDistanceERMVariables(h, x, lambda);
            
            % Calculate the variance of the ERM
            sigma = ones(obj.numContacts, 1) * obj.options.heightVariance;
            dsigma = zeros(obj.numContacts,1 + numel(x) + numel(lambda));
            
            % Scale the distance function
            phi = obj.options.distanceScaling * phi;
            dphi = obj.options.distanceScaling * dphi;
            
            % Calculate the erm Cost
            [f, df] = obj.ermCost(lambda_s, phi, sigma);
            
            % Calculate the total differential cost
            df = df * [dlambda_s; dphi; dsigma];           
            if obj.options.ermSpace == RobustContactImplicitTrajectoryOptimizer.ERM_LOGSPACE
                ep = 1;
                df = 1./(f + ep) .* df;
                f = log(f + ep);
            end
            % Sum over all contacts
            f = obj.options.contactCostMultiplier * sum(f, 1);
            df = obj.options.contactCostMultiplier * sum(df, 1);
            
        end
        function [f, df] = normalDistanceERMSeparated(obj, w)
            %% normalDistanceERMSeparated: Wrapper Function for the ERM Cost when using the Separated Form
            %
            %   Arguments:
            %       OBJ:
            %       z: the slack variables, ordered as [phi, lambdaN]
            %
            %   Return Values
            %       f: ERM Cost evaluation
            %       df: Gradient of ERM Cost wrt z
            
            % ERM Variables
            mu = w(1:obj.numContacts);
            x = w(obj.numContacts+1:end);
            sigma = ones(obj.numContacts,1) * obj.options.heightVariance;  
            % ERM Variable Gradients (wrt z)
            dmu = [eye(numel(mu)), zeros(obj.numContacts, numel(x))];
            dx = [zeros(obj.numContacts, numel(mu)), eye(numel(x))];
            dsigma = zeros(obj.numContacts, numel(w));
            % Calculate the ERM Cost
            [f, df] = obj.ermCost(x, mu, sigma);
            df = df * [dx; dmu; dsigma];
            % Sum over the ERM Cost terms
            f = sum(f, 1);
            df = sum(df, 1);
        end
        function [f, df] = normalDistanceERMSeparatedConstraint(obj,y)
            %% normalDistanceERMSeparatedConstraint: Constraint Enforcing the Slack Variables equal the regular decision variables
            %
            %   normalDistanceERMSeparatedConstraint returns the constraint
            %   evaluation which enforces the ERM Slack variables be equal
            %   to their true values.
            %
            %   Arguments:
            %       y = [h, x, lambda, w]
            
            h = y(1);
            x = y(2:1+obj.plant.getNumStates);
            lambda = y(2+obj.plant.getNumStates:1 + obj.plant.getNumStates + obj.numContacts);
            w = y(2+obj.plant.getNumStates+obj.numContacts:end);
            % Get the ERM Variables
            [phi, lambda_s, dphi, dlambda_s] = obj.getNormalDistanceERMVariables(h, x, lambda);
            % The constraint is the difference between the slack and ERM
            % variables
            f_p = w(1:numel(phi)) - phi;
            f_l = w(numel(phi)+1:end) - lambda_s;
            % Gradients
            df_p = [-dphi, eye(numel(f_p)), zeros(numel(f_p))];
            df_l = [-dlambda_s, zeros(numel(f_l)), eye(numel(f_l))];
            % Combine the results
            f = [f_p; f_l];
            df = [df_p; df_l];
        end
        function [f, df] = normalDistanceERMGradientConstraint(obj, y)
            
            % Get the variables
            h = y(1);
            x = y(2:1+obj.plant.getNumStates);
            lambda = y(2+obj.plant.getNumStates:1 + obj.plant.getNumStates + obj.numContacts);
            % Get the ERM Variables
            [mu, x, dmu, dx, d2mu, d2x] = obj.getNormalDistanceERMVariables(h, x, lambda);
            sigma = ones(obj.numContacts, 1) * obj.options.heightVariance;
            dsigma = zeros(obj.numContacts, numel(y));
            d2sigma = zeros(obj.numContacts, numel(y), numel(y));
            % Get the ERM Gradient and Hessian
            [g, dg, Hg] = obj.ermCost(x, mu, sigma);
            % Collect the Gradient and Hessian into a single vector / matrix
            dz = [dx; dmu; dsigma];
            f = dg * dz;                                    % The gradient is the constraint output
            df = zeros(numel(g), numel(y), numel(y));       % This is the Hessian
            % Combine the parameter hessians
            Hz = zeros(3*numel(g), numel(y), numel(y));
            Hz(1:numel(g),:,:) = d2x;
            Hz(numel(g)+1:2*numel(g),:,:) = d2mu;
            Hz(2*numel(g) + 1:end, :,:) = d2sigma;       
            % Calculate the total ERM Hessian
            for n = 1:numel(g)
               df(n,:,:) = dz' * squeeze(Hg(n,:,:)) * dz + squeeze(sum(dg(n,:)' .* Hz, 1));  
            end
            % Reshape the outputs
            f = reshape(f', [numel(g)*numel(y), 1]);
            df = reshape(permute(df,[2,1,3]), [numel(g)*numel(y), numel(y)]);
        end
        %% ---------------- SLIDING VELOCITY --------------- %%
        function [f, df] = slidingVelocityDefect(obj, x1, lambda)
            %% TODO FIX JACOBIAN DERIVATIVE
            
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
            skip = 2 + obj.numFriction;
            % Sliding velocity slack
            G_idx = skip:skip:nL;
            gamma = lambda(G_idx);
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
            f = f(:);
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
    
    methods (Access = protected)
        function obj = addDistanceERM(obj)
            %% addDistanceERM: add the ERM cost for uncertain terrain height to the problem
            
            grad_level = 1; % Gradients are provided
            switch obj.options.ermMode
                case RobustContactImplicitTrajectoryOptimizer.ERM_OBJECTIVE
                    % Add ERM cost for distance
                    
                    % Create the distance objective
                    distance = FunctionHandleObjective(1 + obj.plant.getNumStates() + obj.numContacts, @obj.normalDistanceERMCost, grad_level);
                    distance = distance.setName(sprintf('DistanceERMCost'));
                    
                    % Add the objective at every point
                    for i = 1:obj.N-1
                        distanceIdx = {obj.h_inds(i); obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds,i)};
                        obj = obj.addCost(distance, distanceIdx);
                    end
                case RobustContactImplicitTrajectoryOptimizer.ERM_SEPARATED
                    % Add ERM cost for distance using slack variables in
                    % the ERM objective
                    distance = FunctionHandleObjective(2*obj.numContacts, @obj.normalDistanceERMSeparated, grad_level);
                    distance = distance.setName(sprintf('DistanceERMCost'));
                    cstr_dim = 1 + obj.plant.getNumStates() + obj.numContacts + 2*obj.numContacts;
                    cstr = FunctionHandleConstraint(zeros(2*obj.numContacts,1), zeros(2*obj.numContacts,1), cstr_dim, @obj.normalDistanceERMSeparatedConstraint);
                    cstr = cstr.setName(sprintf('DistanceERMSlackConstraint'));
                    % Add additional variables to the problem
                    [obj, inds] = obj.addDecisionVariable(2*(obj.N-1)*obj.numContacts); %This adds a lot of decision variables
                    obj.distance_slacks = reshape(inds, 2*obj.numContacts, obj.N-1);
                    % Set the constraints
                    for i = 1:obj.N - 1
                       % Add the distance cost
                       obj = obj.addCost(distance, obj.distance_slacks(:,i));
                       % Add the slack constraint
                       slackIdx = [obj.h_inds(i); obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds, i); obj.distance_slacks(:,i)];
                       obj = obj.addConstraint(cstr, slackIdx);
                    end
                case RobustContactImplicitTrajectoryOptimizer.ERM_CONSTRAINT
                    if obj.options.ermScaling == RobustContactImplicitTrajectoryOptimizer.INERTIA
                       error('Inertia Scaling is currently not supported with the ERM Gradient Constraint method'); 
                    end
                    % Add the gradient of the ERM function as a constraint
                    cstr_dim = 1 + obj.plant.getNumStates + obj.numContacts;
                    distance = FunctionHandleConstraint(zeros(obj.numContacts * cstr_dim, 1), zeros(obj.numContacts * cstr_dim, 1), cstr_dim, @obj.normalDistanceERMGradientConstraint);
                    distance = distance.setName(sprintf('DistanceERMGradientConstraint'));
                    % Add the gradient as a constraint at every point
                    for i = 1:obj.N-1
                       distance_idx = [obj.h_inds(i); obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds,i)];
                       obj = obj.addConstraint(distance, distance_idx);
                    end
                case RobustContactImplicitTrajectoryOptimizer.ERM_COMBINED
                     % Create the distance objective
                    distance = FunctionHandleObjective(1 + obj.plant.getNumStates() + obj.numContacts, @obj.normalDistanceERMCost, grad_level);
                    distance = distance.setName(sprintf('DistanceERMCost'));                   
                    % Add the objective at every point
                    for i = 1:obj.N-1
                        distanceIdx = {obj.h_inds(i); obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds,i)};
                        obj = obj.addCost(distance, distanceIdx);
                    end
                    % add in the complementarity constraint on the normal
                    % distance
                    obj = obj.addDistanceConstraint();
                otherwise
                    error('ERM MODE not recognized');
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
            if obj.options.erm_mode == RobustConactImplictTrajectoryOptimizer.ERM_COMBINED
               obj = obj.addFrictionConstraint(); 
            end
        end
        function obj = addSlidingConstraint(obj)
            %% addSlidingConstraint: adds the complementarity constraint for sliding velocity to the problem
            % Create a Comlementarity Constraint for sliding
            sliding = RelaxedNonlinearComplementarityConstraint(@obj.slidingVelocityConstraint, obj.plant.getNumStates()+2*obj.numContacts, obj.numContacts*obj.numFriction, obj.options.nlcc_mode, obj.options.compl_slack);
            for n = 1:length(sliding.constraints)
                sliding.constraints{n} = sliding.constraints{n}.setName(sprintf('SlidingVelocityConstraint'));
            end
            % Add the constraint at every knot point
            for i = 1:obj.N-1
                % Order the indices such that the argument is 
                %   [x, lambdaN, gamma, lambdaT]
                if obj.options.nlcc_mode == 5
                    slidingIdx = [obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds, i); obj.lambda_inds(obj.gamma_inds, i); obj.lambda_inds(obj.tangent_inds, i); obj.relax_inds(i)];
                else
                    slidingIdx = [obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds, i); obj.lambda_inds(obj.gamma_inds, i); obj.lambda_inds(obj.tangent_inds, i)];
                end
                obj = obj.addConstraint(sliding, slidingIdx);
            end
        end
        function obj = addFrictionConstraint(obj)
            %% addFrictionConstraint: adds the complementarity constraint for the friction cone to the optimization problem
            % Create a complementarity constraint for the friction cone
             friction = RelaxedNonlinearComplementarityConstraint(@obj.frictionConeConstraint, obj.numContacts*(1+obj.numFriction), obj.numContacts, obj.options.nlcc_mode, obj.options.compl_slack);
             for n = 1:length(friction.constraints)
                friction.constraints{n} = friction.constraints{n}.setName(sprintf('FrictionConeConstraint')); 
             end
             % Add the constraint to every knot point                   
             for i = 1:obj.N-1
                 % Reshape all the lambda_inds so the forces are grouped together as [normal, tangential, slack]
                 if obj.options.nlcc_mode == 5
                     frictionIdx =  [obj.lambda_inds(obj.normal_inds,i); obj.lambda_inds(obj.tangent_inds,i); obj.lambda_inds(obj.gamma_inds,i); obj.relax_inds(i)];
                 else
                     frictionIdx =  [obj.lambda_inds(obj.normal_inds,i); obj.lambda_inds(obj.tangent_inds,i); obj.lambda_inds(obj.gamma_inds,i)];
                 end
                 obj = obj.addConstraint(friction, frictionIdx);
             end
        end
        function obj = addDistanceConstraint(obj)
            %% addDistanceConstraint: adds the complementarity constraint for normal distance to the optimization problem 
            % Create a complementarity constraint for the normal distance
            distance = RelaxedNonlinearComplementarityConstraint(@obj.normalDistanceConstraint, obj.plant.getNumStates(), obj.numContacts, obj.options.nlcc_mode, obj.options.compl_slack);
            for n = 1:length(distance.constraints)
                 distance.constraints{n} = distance.constraints{n}.setName(sprintf('NormalDistanceConstraint'));
            end
            % Add the constraint at every knot point
            for i = 1:obj.N-1
                % For distance, we only need [x, lambdaN];
                if obj.options.nlcc_mode == 5
                    distanceIdx = [obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds,i); obj.relax_inds(i)];
                else
                    distanceIdx = [obj.x_inds(:,i+1); obj.lambda_inds(obj.normal_inds,i)];
                end
                obj = obj.addConstraint(distance, distanceIdx);
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
        function [f, df, Hf] = ermCostGaussian(x, mu, sigma)
            
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