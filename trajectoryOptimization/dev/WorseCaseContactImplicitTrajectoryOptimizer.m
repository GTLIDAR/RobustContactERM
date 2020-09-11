classdef WorseCaseContactImplicitTrajectoryOptimizer < ContactImplicitTrajectoryOptimizer
    
   properties
      friction_idx;             % Index for finding the friction coefficient in the decision variables. 
   end
    
    methods
        function obj = WorseCaseContactImplicitTrajectoryOptimizer(plant, N, duration, options)
            
            
            % Create an options structure
            if nargin < 4 || isempty(options)
                options = struct();
            end
            
            % Check the options for friction bounds
            if ~isfield(options, 'lower_friction')
                options.lower_friction = 0.3;       %Lower bound for the friction coefficient
            end
            if ~isfield(options, 'upper_friction')
                options.upper_friction = 0.7;       %Upper bound for the friction coefficient
            end
            
            % Pass construction to the parent class
            obj = obj@ContactImplicitTrajectoryOptimizer(plant, N, duration, options); 
        end
        
        % Overload the Setup method to add the friction coeffiecient as a
        % decision variable
        function obj = setupVariables(obj, N)
           %% SETUP VARIABLES: Overload setupVariables in ContactImplicitTrajectoryOptimizer
           %    to include the friction coefficient as a decision variable   
           
           % Use the parent method first
           obj = setupVariables@ContactImplicitTrajectoryOptimizer(obj, N);
           % Add the friction coefficient as a decision variable
           [obj, obj.friction_idx] = obj.addDecisionVariable(1);
        end
        
        function [xtraj, utraj, ltraj, mu, z, F, info] = solveTraj(obj, t_init, traj_init)
            % Solve the problem using the parent method
            [xtraj, utraj, ltraj, z, F, info] = solveTraj@ContactImplicitTrajectoryOptimizer(obj, t_init, traj_init);
            % Get the value of the friction coefficient
            mu = z(obj.friction_idx);
        end
        
        function obj = addContactConstraints(obj)
            
            % Add the complementarity constraints for contact
            nX = obj.plant.getNumStates();
            nL = obj.numContacts * (2 + obj.numFriction);
            nlc_cstr = NonlinearComplementarityConstraint_original(@obj.contact_constraint_fun, nX + 1, nL, obj.options.nlcc_mode, obj.options.compl_slack);
            worst_cost = FunctionHandleObjective(nL+1, @obj.frictionConeCost, 1);
            worst_cstr = FunctionHandleConstraint(zeros(obj.numContacts,1),inf(obj.numContacts,1),nL, @obj.frictionConeConstraint, 1);
            for i = 1:obj.N - 1
                % Add in the nonlinear complementarity constraints
                
                obj = obj.addConstraint(nlc_cstr, [obj.x_inds(:,i+1); obj.friction_idx; obj.lambda_inds(:,i)]);
                
                % Add in the worst-case objective and constraints on
                % friction coefficient
                obj = obj.addCost(worst_cost, [obj.lambda_inds(:,i); obj.friction_idx]);
                
                % Add in the constraint that the cost be nonnegative
                obj = obj.addConstraint(worst_cstr, obj.lambda_inds(:,i));
                
            end
            
            % Add bounding box constraints on the friction coefficient
            cstr = BoundingBoxConstraint(obj.options.lower_friction, obj.options.upper_friction);
            obj = obj.addBoundingBoxConstraint(cstr, obj.friction_idx);
        end
        
        function z0 = getInitialVars(obj, t_init, traj_init)
            %% getInitialVars: Get the initial decision variables from the trajectories
            %
            %    getInitialVars extends the previous implemention to add
            %    initialization for the friction coefficient
            
            z0 = getInitialVars@ContactImplicitTrajectoryOptimizer(obj, t_init, traj_init);
            % Add friction coefficient - the midpoint of the allowable
            % friction coefficients
            z0(obj.friction_idx) = 1/2*(obj.options.lower_friction + obj.options.upper_friction);
        end
        
        % Overload the contact methods
        function [f, df] = contact_constraint_fun(obj, y)
            %% CONTACT_CONSTRAINT_FUN: Function evaluating the contact nonlinear complementarity function
            %
            % CONTACT_CONSTRIANT_FUN implements the function whose values
            % are complementary to the contact forces. For normal forces
            % lambda_N and frictional forces lambda_T:
            %
            %     lambda_N _|_ phi
            %     lambda_T _|_ gamma + D*dq
            %     gamma    _|_ mu*lambda_N - sum(lambda_T)
            %
            % where phi is the signed distance function, gamma is a slack
            % variable, and D*dq is the relative tangential velocity. The
            % output of CONTACT_NCP_CONSTRAINT_FUN is the direct sum of phi
            % and gamma+D*dq.
            %
            % Input ordering:
            %     The input to contact_ncp_constraint_fun is ordered in a
            %     particular way. The input Y stores the following
            %     information:
            %         Y = [X; MU; LAMBDA; G]
            %     where X is the state of the system, G is the collection of
            %     slack variables gamma, and LAMBDA = [ lambda_N1, lambda_T1,
            %     ..] is the contact forces, alternating between normal and
            %     tangential components.
            
            nQ = obj.plant.getNumPositions();
            nV = obj.plant.getNumVelocities();
            % Get the configuration and velocity
            q = y(1:nQ,:);
            v = y(nQ+1:nQ+nV,:);
            % Get the contact variables
            mu = y(nQ + nV + 1,:);      % Friction coefficient
            z = y(nQ + nV + 2:end, :);
            nZ = numel(z);
            skip = 2 + obj.numFriction;
            T_idx = 1:nZ;
            % Get the normal force
            N_idx = 1:skip:nZ;
            lambda_N = z(N_idx);
            % Get the velocity slack variables
            G_idx = skip:skip:nZ;
            gamma = z(G_idx);
            % Get the tangential force
            T_idx([N_idx,G_idx]) = [];
            lambda_T = z(T_idx);
            
            % Get the contact conditions
            [phi, ~, ~, ~, ~, ~, ~, ~, N, D, ~, dD] = obj.plant.contactConstraints(q, false, obj.options.active_collision_options);
            % Reshape the tangential force basis for easier use
            D = cat(3,D{:});
            D = permute(D,[3,2,1]);
            dD = cat(3,dD{:});
            dD = reshape(dD, [obj.numContacts, nQ, nQ, size(dD, 3)]);
            dD = permute(dD, [4, 2, 3, 1]);
            % Initialize the complementarity function
            f = zeros(nZ, 1);
            df = zeros(nZ, numel(y));
            
            % Distance | Normal complementarity
            f(N_idx,:) = phi;       % Distance function
            df(N_idx,1:nQ) = N;     % Normal vector
                
            for i = 1:obj.numContacts
                % Slack | Cone complementarity
                f(G_idx(i),:) = lambda_N(i)*mu - sum(lambda_T(obj.numFriction*(i-1)+1 : obj.numFriction*i));
                % Derivative wrt normal force
                df(G_idx(i),nQ + nV + 1 + N_idx(i)) = mu;
                % Derivative wrt the friction coefficient
                df(G_idx(i), nQ + nV + 1) = lambda_N(i);
                % Derivative wrt tangential force
                df(G_idx(i),nQ + nV + 1+ T_idx(obj.numFriction*(i-1)+1):nQ + nV + 1 + T_idx(obj.numFriction*i)) = -1;
                % Velocity | Tangential complementarity 
                f(T_idx(obj.numFriction*(i-1)+1):T_idx(obj.numFriction*i),:) = gamma(i) + D(:,:,i)*v;
                % Derivative wrt configuration
                df(T_idx(obj.numFriction*(i-1)+1):T_idx(obj.numFriction*i), 1:nQ) = squeeze(sum(dD(:,:,:,i).*v',2));
                % Derivative wrt velocity
                df(T_idx(obj.numFriction*(i-1)+1):T_idx(obj.numFriction*i), nQ+2:nQ+nV+1) = D(:,:,i);
                % Derivative wrt gamma
                df(T_idx(obj.numFriction*(i-1)+1):T_idx(obj.numFriction*i), nQ+nV+G_idx(i)+1) = 1;
            end
        end
        %% --------------------- WORSE CASE CONSTRAINTS ---------------- %%
        function [f, df] = frictionConeDefect(obj, y)
            % Stored as variables:
            %   y = [lambda, mu]
            
            mu = y(end);
            nZ = numel(y)-1;
            % Pull the contact forces
            skip = 2 + obj.numFriction;
            T_idx = 1:nZ;
            % Get the normal force
            N_idx = 1:skip:nZ;
            lambda_N = y(N_idx);
            % Get the velocity slack variables
            G_idx = skip:skip:nZ;
            % Get the tangential force
            T_idx([N_idx,G_idx]) = [];
            lambda_T = y(T_idx);
            
            % Calculate the friction cone defect
            % Create a matrix to select and sum up the frictional forces
            % for each contact
            S = zeros(obj.numContacts, obj.numContacts * obj.numFriction);
            for n = 1:obj.numContacts
                S(n,obj.numFriction * (n-1) + 1: obj.numFriction * n) = 1;
            end
            % Calculate the friction cone defect.
            f = mu*lambda_N - S * lambda_T;
            % Now calculate the gradients
            df = zeros(numel(f), numel(y));
            % Derivative wrt normal force
            df(:,N_idx) = eye(obj.numContacts) * mu;
            % Derivative wrt tangential force
            df(:,T_idx) = -S;
            % Derivative wrt friction coefficient
            df(:,end) = lambda_N;
        end
        
        function [f, df] = frictionConeCost(obj, y)
           % For the worst case scenario, we negate the friction cone
           % defect to get the worst case friction coefficient (to maximize
           % the friction coefficient, we minimize the negative defect).
           
           % The y-variables are stored as [lambda, mu]
           [f, df] = obj.frictionConeDefect(y);
           f = -sum(f, 1);
           df = -sum(df, 1);
        end
        function [f, df] = frictionConeConstraint(obj, y)
            % We must also ensure that all the friction cone constraints
            % are feasible. So the friction cone with the lowest value must
            % also be feasible. Here, y is simply y = [lambda]
            [f, df] = obj.frictionConeDefect([y(:); obj.options.lower_friction]);
        end
    end
    
    
    
    
end