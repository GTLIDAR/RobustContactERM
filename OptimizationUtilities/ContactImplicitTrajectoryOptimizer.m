classdef ContactImplicitTrajectoryOptimizer < DirectTrajectoryOptimization
   
    properties
        numContacts;    % Number of contacts
        numFriction;    % Number of friction basis vectors per contact
        lambda_inds;    % Indices for finding the contact variables inside the decision variables, ordered [normal; tangential; velocity_slacks];
        force_inds;     % Indices for finding the force variables from inside the contact variables, ordered [normal; tangential];
        normal_inds;    % Row indices for the normal forces    
        tangent_inds;   % Row indices for the tangential forces
        gamma_inds;     % Row indices for the sliding velocity slack variables.
        
        slack_inds;     % Indices for the NCC Slack variables 
        relax_inds = [];% Indices for the NCC relaxation variables         
                
        nJl = 0;        % Number of joint limits
        jl_inds =[];    % Joint limit force indices

          
        lincc_mode = 1; % Mode for linear complementarity constraints (Joint Limits)
        lincc_slack = 0;% Slack variable for linear complementarity constraints (Joint Limits)    
    end  
    properties (Constant)
        % INTEGRATION METHODS
        FORWARD_EULER = 1;
        BACKWARD_EULER = 2;
        MIDPOINT = 3;
        SEMI_IMPLICIT = 4;
        % TIME CONSTRAINT METHODS
        FREETIME = 1;
        PAIRWISETIME = 2;
        BOUNDTIME = 3;
    end
    
    methods 
        function obj = ContactImplicitTrajectoryOptimizer(plant, N, duration, options)
            %% CONTACTIMPLICITTRAJECTORYOPTIMIZER: Construtor Method
            %
            %   
            
            % Create the options structure
            if nargin < 4 
               options = struct(); 
            end
            
            % Check/set some options
            if ~isfield(options,'active_collision_options')
                options.active_collision_options.terrain_only = true;
            end
            if ~isfield(options, 'integration_method')
               options.integration_method = ContactImplicitTrajectoryOptimizer.BACKWARD_EULER; 
            end
            if ~isfield(options,'nlcc_mode')
                options.nlcc_mode = 2;
            end
            if ~isfield(options,'compl_slack')
                options.compl_slack = 0;
            end
            if ~isfield(options, 'time_option')
               options.time_option = 1; 
            end
            if ~isfield(options, 'time_constraints')
                options.time_constraints = 1;
            end
            if ~isfield(options, 'duration')
                if isscalar(duration)
                    options.duration = [duration, duration];
                else
                    options.duration = duration;
                end
            end
            if ~isfield(options, 'relax_cost')
               options.relax_cost = 1; 
            end
            % Construct the object
            obj = obj@DirectTrajectoryOptimization(plant, N, duration, options);    
        end     
        function obj = setupVariables(obj, N)
           %% SETUP VARIABLES: Overload setupVariables in DirectTrajectoryOptimization
           %    to include the forces as decision variables    
           
           % Use the parent method first
           obj = setupVariables@DirectTrajectoryOptimization(obj, N);
           % Then set up additional variables
           % Get the number of contact points and the discretization of
           % friction
           nq = obj.plant.getNumPositions;
           [~, normal, tangential] = obj.plant.contactConstraints(zeros(nq,1));
           obj.numContacts = size(normal, 2);
           obj.numFriction = 2*length(tangential);
           
           % Calculate total number of contact forces
           nContactForces = obj.numContacts * (2 + obj.numFriction);
            % Calculate the indices for getting the forces from the decision
           % variable list
           obj.lambda_inds = reshape(obj.num_vars + (1:N * nContactForces), nContactForces, N);
           
           % Add decision variables for the forces
           obj = obj.addDecisionVariable(N * nContactForces);
                   
           % Set the index variables for the individual forces
           skip = 2 + obj.numFriction;
           % Normal force indices
           obj.normal_inds = 1:skip:nContactForces;
           % Sliding velocity slack indices
           obj.gamma_inds = skip:skip:nContactForces;
           % Tangential force indices
           T_idx = 1:nContactForces;
           T_idx([obj.normal_inds,obj.gamma_inds]) = [];
           obj.tangent_inds = T_idx;
           
           % Set index variables to get just the forces
           obj.force_inds = 1:nContactForces;
           obj.force_inds(obj.gamma_inds) = [];
           
           % Set the relaxation variables
           if obj.options.nlcc_mode == 5
               [obj, inds] = obj.addDecisionVariable(N-1);
               obj.relax_inds = inds;
           end
           
           % Add joint limit forces
           obj.nJl = sum(~isinf([obj.plant.joint_limit_min; obj.plant.joint_limit_max]));
           if obj.nJl > 0
               [obj, inds] = obj.addDecisionVariable(obj.nJl*N);
               obj.jl_inds = reshape(inds,[obj.nJl, N]);
           end
           
           
        end
        function [xtraj, utraj, ftraj,z,F, info, infeasible] = solveTraj(obj, t_init, traj_init)
            % Solve the problem using
            [xtraj, utraj, z, F, info, infeasible] = solveTraj@DirectTrajectoryOptimization(obj, t_init, traj_init);
            
            % Pull the contact forces from the decision variable list
            if obj.numContacts > 0
                t = [0; cumsum(z(obj.h_inds))];
                ftraj = PPTrajectory(foh(t, reshape(z(obj.lambda_inds),[],obj.N)));
            else
                ftraj = [];
            end
        end
        function z0 = getInitialVars(obj, t_init, traj_init)
           %% getInitialVars: Get the initial decision variables from the trajectories
           %
           %    getInitialVars extends the previous implemention to add
           %    contact variables to the decision variable list
           
           z0 = getInitialVars@DirectTrajectoryOptimization(obj, t_init, traj_init);
           % Add contact variables
           if obj.numContacts > 0
               if isfield(traj_init, 'lambda')
                   z0(obj.lambda_inds) = traj_init.lambda.eval(t_init);
               end
           end     
           % Add joint limit variables
           if obj.nJl > 0
              if isfield(traj_init, 'jl')
                    z0(obj.jl_inds) = traj_init.jl.eval(t_init);
              end
           end
        end
        function obj = addDynamicConstraints(obj)
            %% addDynamicConstraints: Add the dynamics as constraints to the problem
            %
            %
            
            % Get some parameters about the system
            nX = obj.plant.getNumStates;
            nU = obj.plant.getNumInputs;
            N = obj.N;
            nL = obj.numContacts * (1 + obj.numFriction);
            % Initialize arrays for holding the constraints
            constraints = cell(N-1, 1);
            dyn_inds = cell(N-1, 1);
            
            % Choose the appropriate dynamic constraint
            switch obj.options.integration_method
                case ContactImplicitTrajectoryOptimizer.FORWARD_EULER
                    n_vars = 2*nX + nU + 1 + nL + obj.nJl;
                    cstr = FunctionHandleConstraint(zeros(nX, 1), zeros(nX, 1), n_vars, @obj.forward_constraint_fun);
                case ContactImplicitTrajectoryOptimizer.BACKWARD_EULER
                    n_vars = 2*nX + nU + 1 + nL + obj.nJl;
                    cstr = FunctionHandleConstraint(zeros(nX, 1), zeros(nX, 1), n_vars, @obj.backward_constraint_fun);
                case ContactImplicitTrajectoryOptimizer.MIDPOINT
                    n_vars = 2*nX + 2*nU + 1 + nL + obj.nJl;
                    cstr = FunctionHandleConstraint(zeros(nX, 1), zeros(nX, 1), n_vars, @obj.midpoint_constraint_fun);
                case ContactImplicitTrajectoryOptimizer.SEMI_IMPLICIT
                    n_vars = 2*nX + nU + 1 + nL + obj.nJl;
                    cstr = FunctionHandleConstraint(zeros(nX, 1), zeros(nX, 1), n_vars, @obj.semiimplicit_constraint_fun);
                otherwise
                    error('Unknown Integration Method');
            end
            cstr = cstr.setName('DynamicConstraints');
            % Get the indices necessary for the constraints (indices to
            % pull the states, controls, and forces from the decision
            % variable list).
            if obj.nJl > 0
                % Add the dynamics
                for i = 1:N-1
                    switch obj.options.integration_method
                        case ContactImplicitTrajectoryOptimizer.FORWARD_EULER
                            dyn_inds{i} = {obj.h_inds(i); obj.x_inds(:,i); obj.x_inds(:,i+1); obj.u_inds(:,i); obj.lambda_inds(obj.force_inds,i); obj.jl_inds(:,i)};
                        case ContactImplicitTrajectoryOptimizer.BACKWARD_EULER
                            dyn_inds{i} = {obj.h_inds(i); obj.x_inds(:,i); obj.x_inds(:,i+1); obj.u_inds(:,i); obj.lambda_inds(obj.force_inds,i); obj.jl_inds(:,i)};
                        case ContactImplicitTrajectoryOptimizer.MIDPOINT
                            dyn_inds{i} = {obj.h_inds(i); obj.x_inds(:,i); obj.x_inds(:,i+1); obj.u_inds(:,i); obj.u_inds(:,i+1); obj.lambda_inds(obj.force_inds,i); obj.jl_inds(:,i)};
                        case ContactImplicitTrajectoryOptimizer.SEMI_IMPLICIT
                            dyn_inds{i} = {obj.h_inds(i); obj.x_inds(:,i); obj.x_inds(:,i+1); obj.u_inds(:,i); obj.lambda_inds(obj.force_inds,i); obj.jl_inds(:,i)};
                        otherwise
                            error('Unknown Integration Method');
                    end
                    constraints{i} = cstr;
                    obj = obj.addConstraint(constraints{i}, dyn_inds{i});
                end
                % Add joint limit constraints
                obj = obj.addJointLimitConstraints();
            else
                % Just add dynamics
                for i = 1:N-1
                    switch obj.options.integration_method
                        case ContactImplicitTrajectoryOptimizer.FORWARD_EULER
                            dyn_inds{i} = {obj.h_inds(i); obj.x_inds(:,i); obj.x_inds(:,i+1); obj.u_inds(:,i); obj.lambda_inds(obj.force_inds,i)};
                        case ContactImplicitTrajectoryOptimizer.BACKWARD_EULER
                            dyn_inds{i} = {obj.h_inds(i); obj.x_inds(:,i); obj.x_inds(:,i+1); obj.u_inds(:,i); obj.lambda_inds(obj.force_inds,i)};
                        case ContactImplicitTrajectoryOptimizer.MIDPOINT
                            dyn_inds{i} = {obj.h_inds(i); obj.x_inds(:,i); obj.x_inds(:,i+1); obj.u_inds(:,i); obj.u_inds(:,i+1); obj.lambda_inds(obj.force_inds,i)};
                        case ContactImplicitTrajectoryOptimizer.SEMI_IMPLICIT
                            dyn_inds{i} = {obj.h_inds(i); obj.x_inds(:,i); obj.x_inds(:,i+1); obj.u_inds(:,i); obj.lambda_inds(obj.force_inds,i)};
                        otherwise
                            error('Unknown Integration Method');
                    end
                    constraints{i} = cstr;
                    obj = obj.addConstraint(constraints{i}, dyn_inds{i});
                end
            end
            % Add in the complementarity constriants
            obj = obj.addContactConstraints();
            % Add constraints on the time steps
            if obj.options.time_option == 2
               obj = obj.addTimeStepConstraints(); 
            end
            
        end   
        function obj = addContactConstraints(obj)
            
            nX = obj.plant.getNumStates();
            % NC Constraint with no slack variables - original constraint
            nlc_cstr = RelaxedNonlinearComplementarityConstraint(@obj.contact_constraint_fun, nX, obj.numContacts * (2 + obj.numFriction), obj.options.nlcc_mode, obj.options.compl_slack);
            nlc_cstr.constraints{1} = nlc_cstr.constraints{1}.setName('ContactForceNonneg');
            nlc_cstr.constraints{2} = nlc_cstr.constraints{2}.setName('ContactFuncNonneg');
            nlc_cstr.constraints{3} = nlc_cstr.constraints{3}.setName('ContactCompl');
            if obj.options.nlcc_mode == 5
                % Relaxed NC Constraints
                for i = 1:obj.N - 1
                    obj = obj.addConstraint(nlc_cstr, [obj.x_inds(:,i+1); obj.lambda_inds(:,i); obj.relax_inds(i)]);
                end
                % Add the cost for the relaxation
                cost = FunctionHandleObjective(obj.N-1, @(x) obj.relaxed_nc_cost(x));
                cost = cost.setName('RelaxedNLCCost');
                obj = obj.addCost(cost, obj.relax_inds);
            else
                % Strict NC Constraints
                for i = 1:obj.N - 1
                    obj = obj.addConstraint(nlc_cstr, [obj.x_inds(:,i+1); obj.lambda_inds(:,i)]);
                end
            end
            
        end
        function obj = addJointLimitConstraints(obj)
            %% Add joint limit constraints to the problem
            
            nQ = obj.plant.getNumPositions();
            % Create a linear complementarity constraint for the joint
            % limits
            W = zeros(obj.nJl);
            M = [eye(nQ); -eye(nQ)];
            b = [-obj.plant.joint_limit_min; obj.plant.joint_limit_max];
            inflimit = isinf(b);
            M(inflimit,:) = [];
            b(inflimit) = [];
            cstr = LinearComplementarityConstraint_original(W, b, M, obj.lincc_mode, obj.lincc_slack);
            cstr.constraints{1} = cstr.constraints{1}.setName('JointForceNonneg');
            cstr.constraints{2} = cstr.constraints{2}.setName('JointLimitNonneg');
            cstr.constraints{3} = cstr.constraints{3}.setName('JointLimitCompl');
            for n = 1:obj.N-1
                obj = obj.addConstraint(cstr,[obj.x_inds(1:nQ,n+1); obj.jl_inds(:,n)]);
            end
            
            
        end
        function obj = addRunningCost(obj, running_cost_function, name)
           %% addRunningCost: add the running cost function as the objective
           %
           %    
            nX = obj.plant.getNumStates();
            nU = obj.plant.getNumInputs();
            for i=1:obj.N-1
                switch obj.options.integration_method
                    case ContactImplicitTrajectoryOptimizer.FORWARD_EULER
                        running_cost = FunctionHandleObjective(1+nX+nU, @(h, x, u)obj.reimann_sum(running_cost_function, h, x, u));
                        inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.u_inds(:,i)};
                    case ContactImplicitTrajectoryOptimizer.BACKWARD_EULER
                        running_cost = FunctionHandleObjective(1+nX+nU, @(h, x, u)obj.reimann_sum(running_cost_function, h, x, u));
                        inds_i = {obj.h_inds(i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
                    case ContactImplicitTrajectoryOptimizer.MIDPOINT
                        running_cost = FunctionHandleObjective(1+2*nX+2*nU,...
                            @(h,x0,x1,u0,u1) obj.midpoint_running_fun(running_cost_function,h,x0,x1,u0,u1));
                        inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1)};
                    case ContactImplicitTrajectoryOptimizer.SEMI_IMPLICIT
                        running_cost = FunctionHandleObjective(1+nX+nU, @(h, x, u)obj.reimann_sum(running_cost_function, h, x, u));
                        inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.u_inds(:,i)};
                    otherwise
                        error('Unknown integration method');
                end
                if nargin == 3
                    running_cost = running_cost.setName(name);
                else
                    running_cost = running_cost.setName(sprintf('RunningCost'));
                end
                obj = obj.addCost(running_cost,inds_i);
            end
        end
        function obj = addForceCost(obj, cost_func, name)
            %% addForceCost: Add a cost term on the contact forces
            %
            % addForceCosts allows the user to add a cost objective to
            % the contact forces at every node point in the problem.
            % The cost must have the following syntax:
            %
            %   [g, dg] = force_cost(h, l)
            % 
            % where:
            %   h is the scalar timestep of the current node
            %   l is a vector of contact forces, stored as [normal,
            %   tangential]
            %   g is the scalar cost 
            %   dg is the gradient wrt the timestep h and the forces l
            
            % Total number of contact forces
            nL = obj.numContacts*(1 + obj.numFriction);
            
            % Create and name the cost handle
            cost = FunctionHandleObjective(1 + nL, @(h, l)cost_func(h, l));
            if nargin == 3
                cost = cost.setName(name);
            else
                cost = cost.setName(sprintf('ForceCost'));
            end
            % Add the cost at every node point
            for n = 1:obj.N-1
               inds = {obj.h_inds(n); obj.lambda_inds(obj.force_inds, n)};
               obj = obj.addCost(cost, inds);
            end
        end
        function obj = addDifferenceCost(obj, cost_func, name)
            %% addDifferenceCost: adds a cost on the difference in values between node points
            %
            %   addDifferenceCost adds a cost on the difference between
            %   node points. The differenced cost can be the difference on
            %   controls, states, or forces. The cost must have the syntax:
            %
            %   [g,dg] = cost(h, dx, du, dl)
            %  
            %   where
            %       h is the different in time points between nodes
            %       dx is the difference in states between nodes
            %       du is the difference in controls between nodes
            %       dl is the difference in contact forces between nodes
                        
            % State, control, and force dimensions
            nX = obj.plant.getNumStates();
            nU = obj.plant.getNumInputs();
            nL = obj.numContacts * (1 + obj.numFriction);
            
            % Create a cost handle
            running_cost = FunctionHandleObjective(1 + 2*nX + 2*nU + 2*nL, @(h, x0, x1, u0, u1, l0, l1)obj.differenceCost(cost_func, h, x0, x1, u0, u1, l0, l1));
            if nargin ==3
                running_cost = running_cost.setName(name);
            else
                running_cost = running_cost.setName(sprintf('DifferencedCost'));
            end
            % Add the cost at every node point
            for n = 1:obj.N-1
                inds = {obj.h_inds(n); obj.x_inds(:,n); obj.x_inds(:,n+1); obj.u_inds(:, n); obj.u_inds(:,n+1);...
                    obj.lambda_inds(obj.force_inds,n); obj.lambda_inds(obj.force_inds, n+1)};
                obj = obj.addCost(running_cost, inds);
            end
            
        end
        function obj = addTimeStepConstraints(obj)
            %% Add Timestep Constriants: Adds constraints to the values of the timesteps
            %   addTimeStepConstraints can add the following constraints:
            %       1. No Constriants: each timestep can be any nonnegative
            %       number
            %       2. Pairwise constriants: every two timesteps must have
            %       the same total duration, although that duration is not
            %       specified
            %       3. Bounded Constriants: every timestep must be between
            %       two values calculated from the number of knot points
            %       and the total duration of the problem.
            
            switch obj.options.time_constraints
                case ContactImplicitTrajectoryOptimizer.FREETIME
                    % There are no constraints on time.
                case ContactImplicitTrajectoryOptimizer.PAIRWISETIME
                    % Pairwise time constraints specify that every two
                    % steps must have the same collective duration
                    
                    % Add pairwise constraints over the timesteps
                    Ncol = obj.N-1;
                    Nrow = floor((obj.N-3)/2);
                    % Create a linear constraint matrix
                    A = zeros(Nrow, Ncol);
                    for n = 1:Nrow
                        A(n, 2*(n-1)+1:2*(n+1)) = [1, 1, -1, -1];
                    end
                    if rem(Ncol, 2) ~= 0
                        % Add an additional row if N is odd
                        A = [A; zeros(1, Ncol)];
                        A(end, end - 2:end) = [1, 1, -1];
                    end
                    % Create a linear constraint
                    cstr = LinearConstraint(zeros(1,size(A,1)), zeros(1,size(A,1)), A);
                    cstr = cstr.setName(sprintf('PairwiseTimeConstraints'));
                    % Add the constraint to the problem
                    obj = obj.addConstraint(cstr, obj.h_inds);
                case ContactImplicitTrajectoryOptimizer.BOUNDTIME
                    % Calculate the necessary timestep for the minimum
                    % duration
                    dt_min = obj.options.duration(1)/(obj.N - 1);
                    % Calculate the necessary timestep for the maximum duration
                    dt_max = obj.options.duration(2)/(obj.N - 2);
                    % Let the bounds be 50% the min and 150% the max
                    cstr = BoundingBoxConstraint(0.5*dt_min*ones(1,obj.N-1), 1.5*dt_max*ones(1,obj.N-1));
                    cstr = cstr.setName(sprintf('BoundedTimeConstriants'));
                    % Place the constraint over the timesteps
                    obj = obj.addConstraint(cstr, obj.h_inds);
                otherwise
                    % No Constraints
                    error('Unrecognized timestep constraint method');
            end
        end
        function obj = addFinalCost(obj, cost_func, name)
            %% addFinalCost: Adds a terminal cost to the problem
            %
            %   addFinalCost overloads the method in
            %   DirectTrajectoryOptimization to allow the user to name the
            %   final cost for display. The cost function included must
            %   have the following syntax:
            %
            %   [f, df] = final_cost(h, x)
            %
            %   where 
            %       h is a vector of all timesteps, 
            %       x is the final state, 
            %       f is the value of the final cost
            %       df is the gradient of the final cost with respect to h
            %       and x
            
            % Dimensions
            nX = obj.plant.getNumStates();
            nH = obj.N - 1;
            % Create the handle obejctive
            cost = FunctionHandleObjective(nH + nX, @(h, x)cost_func(h, x));
            % Label the terminal cost
            if nargin == 3
                cost = cost.setName(name);
            else
                cost = cost.setName('TerimalCost');
            end
            % Add the cost to the problem
            obj = obj.addCost(cost, {obj.h_inds; obj.x_inds(:,end)});
        end
    end
    methods 
        function [f, df] = forward_constraint_fun(obj, h, x0, x1, u, lambda, jlambda)
           %% FORWARD_CONSTRAINT_FUN: Enforces dynamics using a Forward-Euler integration
           %
           %    Forward Constraint Function enforces the dynamics
           %    constraint using a forward Euler approximation to the
           %    dynamics. For the differential equation:
           %        dx = g(x)
           %    The forward Euler approximation calculates the value of x
           %    at step k+1 from the value at step k according to:
           %        x[k+1] = x[k] + hg(x[k])
           %
           %    where h is an appropriate step size.
           %
           %    Arguments:
           %        OBJ:    A ContactImplicitTrajectoryOptimizer object
           %        h:      Scalar double, the step size
           %        x0:     nXx1 double, initial state (state at step k)
           %        x1:     nXx1 double, final state (state at step k+1)
           %        u:      nUx1 double, the controls at step k
           %        lambda: nLx1 double, the contact forces at step k
           %
           %    Return values:
           %        f:      nX x 1 double, the dynamic defects
           %        df:     nX x (2*nX+nU+nL+1) double, the derivatives of
           %                the defects with respect to each of the inputs 
           %                (in order)
           %
           % The defects f are calculated according to the forward Euler
           % rule. Specifically:
           %        f = x[k+1] - x[k] - hg(x[k])
           
           % Get some parameters about the system
           nQ = obj.plant.getNumPositions;
           nV = obj.plant.getNumVelocities;
           nU = obj.plant.getNumInputs;
           nL = length(lambda);
           
           % Ensure the velocities and positions have equal number of
           % elements
           assert(nQ == nV);
            
           % Get the configuration and velocity at the endpoints of the
           % interval
           q0 = x0(1:nQ);
           v0 = x0(nQ+1:nQ+nV);
           q1 = x1(1:nQ);
           v1 = x1(nQ+1:nQ+nV);
           
           % FORWARD INTEGRATION: Get the values of the dynamics at the
           % beginning of the interval
           [H, C, B, dH, dC, dB] = obj.plant.manipulatorDynamics(q0, v0);
           % Get the contact jacobian
           [J, dJ]  = obj.evalContactJacobian(q1);
           % Reshape the gradients back into tensor form
           dH = reshape(dH',[nQ, nQ, nQ]);
           dB = reshape(dB',[nQ, nU, nQ]);
           dJ = permute(dJ, [2,1,3]);
           % Calculate the forward dynamics defects (position equation)
           fq = q1 - q0 - h * v0;
           dfq = [-v0,-eye(nQ),-h*eye(nV),eye(nQ), zeros(nQ, nV + nU + nL)];  
           % Calculate the forward dynamics defects (velocity equation)
           fv = H*(v1 - v0) - h*(B*u - C + J'*lambda);
           % Calculate the derivative
           dHv = squeeze(sum(dH .* (v1 - v0)', 2));
           dBu = squeeze(sum(dB .* u', 2));
           dJl = squeeze(sum(dJ .* lambda', 2));
           dfv = [-(B*u - C + J'*lambda), dHv - h*(dBu - dC(:,1:nQ)), -H + h*dC(:,nQ+1:nQ+nV), -h*dJl, H, -h*B, -h*J'];
           % Add joint limits
           if obj.nJl > 0
               Jl = [eye(nQ); -eye(nQ)];
               b = [-obj.plant.joint_limit_min; obj.plant.joint_limit_max];
               Jl(isinf(b),:) = [];
               % Add to the dynamics
               fv = fv - h*Jl'*jlambda;
               % Add to the gradients
               dfv(:,1) = dfv(:,1) - Jl'*jlambda;
               dfv = [dfv, -h*Jl'];
               dfq = [dfq, zeros(size(Jl'))];
           end
           % Combine the defects and the derivatives
           f = [fq; fv];
           df = [dfq; dfv];
        end
        function [f, df] = backward_constraint_fun(obj, h, x0, x1, u, lambda, jlambda)
            %% BACKWARD_CONSTRAINT_FUN: Enforces dynamics using a Backward-Euler integration
            %
            %    Backward Constraint Function enforces the dynamics
            %    constraint using a backward Euler approximation to the
            %    dynamics. For the differential equation:
            %        dx = g(x)
            %    The backward Euler approximation calculates the value of x
            %    at step k+1 from the value at step k according to:
            %        x[k+1] = x[k] + hg(x[k+1])
            %
            %    where h is an appropriate step size.
            %
            %    Arguments:
            %        OBJ:    A ContactImplicitTrajectoryOptimizer object
            %        h:      Scalar double, the step size
            %        x0:     nXx1 double, initial state (state at step k)
            %        x1:     nXx1 double, final state (state at step k+1)
            %        u:      nUx1 double, the controls at step k
            %        lambda: nLx1 double, the contact forces at step k
            %
            %    Return values:
            %        f:      nX x 1 double, the dynamic defects
            %        df:     nX x (2*nX+nU+nL+1) double, the derivatives of
            %                the defects with respect to each of the inputs
            %                (in order)
            %
            % The defects f are calculated according to the backward Euler
            % rule. Specifically:
            %        f = x[k+1] - x[k] - hg(x[k+1])
            
            % Get some parameters about the system
            nQ = obj.plant.getNumPositions;
            nV = obj.plant.getNumVelocities;
            nU = obj.plant.getNumInputs;
            nL = length(lambda);
            
            % Ensure the velocities and positions have equal number of
            % elements
            assert(nQ == nV);
            
            % Get the configuration and velocity at the endpoints of the
            % interval
            q0 = x0(1:nQ);
            v0 = x0(nQ+1:nQ+nV);
            q1 = x1(1:nQ);
            v1 = x1(nQ+1:nQ+nV);
            
            % BACKWARD INTEGRATION: Get the values of the dynamics at the
            % end of the interval
            [H, C, B, dH, dC, dB] = obj.plant.manipulatorDynamics(q1, v1);
            % Get the contact jacobian (always at the end of the interval)
            [J, dJ]  = obj.evalContactJacobian(q1);
            % Reshape the gradients back into tensor form
            dH = reshape(dH',[nQ, nQ, nQ]);
            dB = reshape(dB',[nQ, nU, nQ]);
            dJ = permute(dJ, [2,1,3]);
            % Calculate the forward dynamics defects (position equation)
            fq = q1 - q0 - h * v1;
            dfq = [-v1,-eye(nQ),zeros(nQ,nV),eye(nQ),-h*eye(nV),zeros(nQ,nU + nL)];
            % Calculate the forward dynamics defects (velocity equation)
            fv = H*(v1 - v0) - h*(B*u - C + J'*lambda);
%             fv = H*(v1 - v0) - h*(B*u - C) - J'*lambda;
            % Calculate the derivative
            dHv = squeeze(sum(dH .* (v1 - v0)', 2));
            dBu = squeeze(sum(dB .* u', 2));
            dJl = squeeze(sum(dJ .* lambda', 2));
            dfv = [-(B*u - C + J'*lambda), zeros(nV, nQ), -H, dHv - h*(dBu - dC(:,1:nQ) + dJl), H + h*dC(:,nQ+1:nQ+nV), -h*B, -h*J'];
%             dfv = [-B*u + C, zeros(nV, nQ), -H, dHv - h*(dBu - dC(:,1:nQ) + dJl), H + h*dC(:,nQ+1:nQ+nV), -h*B, -J'];
            % Add joint limits
            if obj.nJl > 0
                Jl = [eye(nQ); -eye(nQ)];
                b = [-obj.plant.joint_limit_min; obj.plant.joint_limit_max];
                Jl(isinf(b),:) = [];
                % Add to the dynamics
                fv = fv - h*Jl'*jlambda;
%                 fv = fv - Jl'*jlambda;
                % Add to the gradients
                dfv(:,1) = dfv(:,1) - Jl'*jlambda;
                dfv = [dfv, -h*Jl'];
%                 dfv = [dfv, -Jl'];
                dfq = [dfq, zeros(size(Jl'))];
            end
            % Combine the defects and the derivatives
            f = [fq; fv];
            df = [dfq; dfv];
        end
        function [f, df] = midpoint_constraint_fun(obj, h, x0, x1, u0, u1, lambda, jlambda)
            %% MIDPOINT_CONSTRAINT_FUN: Enforces dynamics using Midpoint integration
            %
            %    Midpoint Constraint Function enforces the dynamics
            %    constraint using a midpoint approximation to the
            %    dynamics. For the differential equation:
            %        dx = g(x)
            %    The midpoint calculates the value of x at step k+1 from 
            %    the value at step k according to:
            %        x[k+1] = x[k] + hg(1/2 *(x[k] + x[k+1]))
            %
            %    where h is an appropriate step size.
            %
            %    Arguments:
            %        OBJ:    A ContactImplicitTrajectoryOptimizer object
            %        h:      Scalar double, the step size
            %        x0:     nXx1 double, initial state (state at step k)
            %        x1:     nXx1 double, final state (state at step k+1)
            %        u0:     nUx1 double, the controls at step k
            %        u1:     nUx1 double, the controls at step k+1
            %        lambda: nLx1 double, the contact forces at step k
            %
            %    Return values:
            %        f:      nX x 1 double, the dynamic defects
            %        df:     nX x (2*nX+2*nU+nL+1) double, the derivatives of
            %                the defects with respect to each of the inputs
            %                (in order)
            %
            % The defects f are calculated according to the midpoint rule.
            % Specifically:
            %        f = x[k+1] - x[k] - hg(1/2 * (x[k] + x[k+1]))

            % Get some parameters about the system
            nQ = obj.plant.getNumPositions;
            nV = obj.plant.getNumVelocities;
            nU = obj.plant.getNumInputs;
            nL = length(lambda);
            
            % Ensure the velocities and positions have equal number of
            % elements
            assert(nQ == nV);
            
            % Get the configuration and velocity at the endpoints of the
            % interval
            q0 = x0(1:nQ);
            v0 = x0(nQ+1:nQ+nV);
            q1 = x1(1:nQ);
            v1 = x1(nQ+1:nQ+nV);
            % Calculate the midpoint control
            u = 0.5*(u0 + u1);
            
            % MIDPOINT INTEGRATION: Evaluate the constraints at the
            % midpoint of the interval
            [H, C, B, dH, dC, dB] = obj.plant.manipulatorDynamics(0.5*(q0+q1),0.5*(v0+v1));
            % Get the contact jacobian (always at the end of the interval)
            [J, dJ]  = obj.evalContactJacobian(q1);
            % Reshape the gradients back into tensor form
            dH = reshape(dH',[nQ, nQ, nQ]);
            dB = reshape(dB',[nQ, nU, nQ]);
            dJ = permute(dJ, [2,1,3]);
            % Calculate the forward dynamics defects (position equation)
            fq = q1 - q0 - h * 0.5 * (v0+v1);
            dfq = [-0.5*(v0+v1),-eye(nQ),-h*0.5*eye(nV),eye(nQ),-h*0.5*eye(nV) zeros(nQ, 2*nU + nL)];
            % Calculate the forward dynamics defects (velocity equation)
            fv = H*(v1 - v0) - h*(B*u - C + J'*lambda);
            % Calculate the derivative
            dHv = squeeze(sum(dH .* (v1 - v0)', 2));
            dBu = squeeze(sum(dB .* u', 2));
            dJl = squeeze(sum(dJ .* lambda', 2));
            dfv = [-(B*u - C + J'*lambda), 0.5*(dHv - h*(dBu - dC(:,1:nQ))), -H + h*0.5*dC(:,nQ+1:nQ+nV),...
                0.5*(dHv - h*(dBu - dC(:,1:nQ))) - h*dJl, H + h*0.5*dC(:,nQ+1:nQ+nV), -h*0.5*B,-h*0.5*B, -h*J'];
            % Add joint limits
            if obj.nJl > 0
                Jl = [eye(nQ); -eye(nQ)];
                b = [-obj.plant.joint_limit_min; obj.plant.joint_limit_max];
                Jl(isinf(b),:) = [];
                % Add to the dynamics
                fv = fv - h*Jl'*jlambda;
                % Add to the gradients
                dfv(:,1) = dfv(:,1) - Jl'*jlambda;
                dfv = [dfv, -h*Jl'];
                dfq = [dfq, zeros(size(Jl'))];
            end
            % Combine the defects and the derivatives
            f = [fq; fv];
            df = [dfq; dfv];
        end
        function [f, df] = semiimplicit_constraint_fun(obj, h, x0, x1, u, lambda, jlambda)
            %% SEMIIMPLICIT_CONSTRAINT_FUN: Enforces dynamics using Semi-Implicit Euler integration
            %
            %    Semi-Implicit Constraint Function enforces the dynamics
            %    constraint using a semi-implicit Euler approximation to the
            %    dynamics. For the set of differential equations
            %        dx = g(v)
            %        dv = f(x)
            %    The semi-implicit Euler approximation calculates the 
            %    values of x and v at step k+1 from the value at step k 
            %    according to:
            %        v[k+1] = v[k] + hf(x[k])
            %        x[k+1] = x[k] + hg(v[k+1])
            %
            %    where h is an appropriate step size.
            %
            %    Arguments:
            %        OBJ:    A ContactImplicitTrajectoryOptimizer object
            %        h:      Scalar double, the step size
            %        x0:     nXx1 double, initial state (state at step k)
            %        x1:     nXx1 double, final state (state at step k+1)
            %        u:      nUx1 double, the controls at step k
            %        lambda: nLx1 double, the contact forces at step k
            %
            %    Return values:
            %        f:      nX x 1 double, the dynamic defects
            %        df:     nX x (2*nX+nU+nL+1) double, the derivatives of
            %                the defects with respect to each of the inputs
            %                (in order)
            
            % Get some parameters about the system
            nQ = obj.plant.getNumPositions();
            nV = obj.plant.getNumVelocities();
            nU = obj.plant.getNumInputs();
            nL = length(lambda);
            
            % Ensure the velocities and positions have equal number of
            % elements
            assert(nQ == nV);
            
            % Get the configuration and velocity at the endpoints of the
            % interval
            q0 = x0(1:nQ);
            v0 = x0(nQ+1:nQ+nV);
            q1 = x1(1:nQ);
            v1 = x1(nQ+1:nQ+nV);
            
            % SEMI-IMPLICIT INTEGRATION: Get the values of the dynamics at the
            % beginning of the interval
            [H, C, B, dH, dC, dB] = obj.plant.manipulatorDynamics(q0, v0);
            % Get the contact jacobian
            [J, dJ]  = obj.evalContactJacobian(q1);
            % Reshape the gradients back into tensor form
            dH = reshape(dH',[nQ, nQ, nQ]);
            dB = reshape(dB',[nQ, nU, nQ]);
            dJ = permute(dJ, [2,1,3]);
            % Calculate the backward dynamics defects (position equation)
            fq = q1 - q0 - h * v1;
            dfq = [-v0,-eye(nQ),zeros(nQ,nV) ,eye(nQ), -h*eye(nV), zeros(nQ, nU + nL)];
            % Calculate the forward dynamics defects (velocity equation)
            fv = H*(v1 - v0) - h*(B*u - C + J'*lambda);
            % Calculate the derivative
            dHv = squeeze(sum(dH .* (v1 - v0)', 2));
            dBu = squeeze(sum(dB .* u', 2));
            dJl = squeeze(sum(dJ .* lambda', 2));
            dfv = [-(B*u - C + J'*lambda), dHv - h*(dBu - dC(:,1:nQ)), -H + h*dC(:,nQ+1:nQ+nV), -h*dJl, H, -h*B, -h*J'];
            % Add joint limits
            if obj.nJl > 0
                Jl = [eye(nQ); -eye(nQ)];
                b = [-obj.plant.joint_limit_min; obj.plant.joint_limit_max];
                Jl(isinf(b),:) = [];
                % Add to the dynamics
                fv = fv - h*Jl'*jlambda;
                % Add to the gradients
                dfv(:,1) =dfv(:,1) - Jl'*jlambda;
                dfv = [dfv, -h*Jl'];
                dfq = [dfq, zeros(size(Jl'))];
            end
            % Combine the defects and the derivatives
            f = [fq; fv];
            df = [dfq; dfv];
        end
        function [J, dJ] = evalContactJacobian(obj,q)
            
            % Get the normal and tangential components of the Jacobian
            [~,~, ~, ~, ~, ~, ~, ~, N, D, dN, dD] = obj.plant.contactConstraints(q, false, obj.options.active_collision_options);
            % Concatenate the components of the Jacobian together
            % Note: we must concatenate these in the correct order
            J = cell(1,obj.numContacts);
            dJ = cell(1,obj.numContacts);
            % Get rid of the cell arrays to facilitate concatenation
            % Note: The friction vectors are first collected across the
            % different contact points, and then collected across the basis
            % vectors. That is, if there are 3 contact points and 2
            % friction bases, then rows 1-3 are the positive friction
            % vectors and rows 4-6 are the negative friction vectors.
            D = cat(1,D{:});            
            % Reshape the gradients
            dN = reshape(dN',[numel(q), numel(q), obj.numContacts]);
            dN = permute(dN,[3,1,2]);
            for n = 1:length(dD)
                dD{n} = reshape(dD{n}',numel(q), numel(q), obj.numContacts);
                dD{n} = permute(dD{n},[3,1,2]);
            end
            dD = cat(1,dD{:});
            %% TODO: FIX THE JACOBIAN DERIVATIVES FOR MULTI-CONTACT
            for k = 1:obj.numContacts
                % Jacobian for the kth contact
                J{k} = [N(k,:); D(k:obj.numContacts:end,:)];
                % Jacobian derivative for the kth contact
                dJ{k} = cat(1,dN(k,:,:),dD(k:obj.numContacts:end,:,:));
            end
            J = cat(1, J{:});
            dJ = cat(1,dJ{:});
        end
        function [g, dg] = midpoint_running_fun(obj, cost_fun, h, x0, x1, u0, u1)
            
            nX = obj.plant.getNumStates();
            nU = obj.plant.getNumInputs();
            [g, dh] = cost_fun(h, 0.5*(x0 + x1), 0.5*(u0 + u1));
            % Include the timestep in the integration
            g = h*g;
            dh = h * dh;
            
            dg = [dh(:,1) + g, 0.5*dh(:,2:1 + nX), 0.5*dh(:, 2:1 + nX), 0.5 * dh(:, nX+2 : nX+nU+1), 0.5 * dh(:, nX+2: nX + nU + 1)];
        end
        function [g, dg] = reimann_sum(obj, cost_fun, h, x0, u0)
            
            [f, df] = cost_fun(h, x0, u0);
            g = h * f;
            dg = h * df;
            dg(:,1) = dg(:,1) + f;
        end
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
            %         Y = [X; LAMBDA; G]
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
            z = y(nQ + nV + 1:end, :);
            nZ = numel(z);
            lambda_N = z(obj.normal_inds);
            gamma = z(obj.gamma_inds);
            lambda_T = z(obj.tangent_inds);
            
            % Get the contact conditions
            [phi, ~, ~, ~, ~, ~, ~, mu, N, D, ~, dD] = obj.plant.contactConstraints(q, false, obj.options.active_collision_options);
            % Reshape the tangential force basis for easier use
            D = cat(1,D{:});
            for n = 1:length(dD)
                dD{n} = reshape(dD{n}', [nQ, nQ, obj.numContacts]);
                dD{n} = permute(dD{n}, [3, 1, 2]);
            end
            dD = cat(1,dD{:});
            % Initialize the complementarity function
            f = zeros(nZ, 1);
            df = zeros(nZ, nQ + nV + obj.numContacts*(2 + obj.numFriction));
            
            % Distance | Normal complementarity
            f(obj.normal_inds,:) = phi;       % Distance function
            df(obj.normal_inds,1:nQ) = N;     % Normal vector
                
            for i = 1:obj.numContacts
                % Slack | Cone complementarity
                f(obj.gamma_inds(i),:) = lambda_N(i)*mu(i) - sum(lambda_T(obj.numFriction*(i-1)+1 : obj.numFriction*i));
                % Derivative wrt normal force
                df(obj.gamma_inds(i),nQ + nV + obj.normal_inds(i)) = mu(i);
                % Derivative wrt tangential force
                df(obj.gamma_inds(i),nQ + nV + obj.tangent_inds(obj.numFriction*(i-1)+1):nQ + nV + obj.tangent_inds(obj.numFriction*i)) = -1;
                % Velocity | Tangential complementarity 
                f(obj.tangent_inds(obj.numFriction*(i-1)+1):obj.tangent_inds(obj.numFriction*i),:) = gamma(i) + D(i:obj.numContacts:end,:)*v;
                % Derivative wrt configuration
                df(obj.tangent_inds(obj.numFriction*(i-1)+1):obj.tangent_inds(obj.numFriction*i), 1:nQ) = squeeze(sum(dD(i:obj.numContacts:end,:,:).*v',2));
                % Derivative wrt velocity
                df(obj.tangent_inds(obj.numFriction*(i-1)+1):obj.tangent_inds(obj.numFriction*i), nQ+1:nQ+nV) = D(i:obj.numContacts:end,:);
                % Derivative wrt gamma
                df(obj.tangent_inds(obj.numFriction*(i-1)+1):obj.tangent_inds(obj.numFriction*i), nQ+nV+obj.gamma_inds(i)) = 1;
            end
        end
        function [f, df] = differenceCost(obj, cost_fun, h, x0, x1, u0, u1, l0, l1)
            %% DifferencedCost: Private function for handling costs on differenced terms
            %
            %   differencedCosst implements a general handle for including
            %   costs on the difference between successive nodes in the
            %   optimization. differencedCost adds a cost on the difference
            %   between successive states, controls, and contact forces.
            
            % Get the number of states and inputs
            nX = obj.plant.getNumStates();
            nU = obj.plant.getNumInputs();
            % Take the difference between successive node points
            dx = x1 - x0;
            du = u1 - u0;
            dl = l1 - l0;
            % Calculate the cost on the differences
            [f,dg] = cost_fun(h, dx, du, dl);
            % Reshape the derivaties (apply chain rule)
            df = [dg(1), -dg(2:1+nX), dg(2:1+nX), -dg(nX+2:nX+1+nU), dg(nX+2:nX+1+nU), ...
               -dg(nX + nU + 2:end), dg(nX + nU + 2:end)];
        end
    end
    methods
        function obj = enableCostDisplay(obj)
            %% ENABLECOSTDISPLAY: Display Each of the Costs and Constraints associated with the problem
            
            % Clear and persistent variables
           obj.printCostsAndConstraints(); 
           % Add the display function to the list of displays
           obj = obj.addDisplayFunction(@(z)obj.printCostsAndConstraints(z));
        end      
        function obj = printCostsAndConstraints(obj, z)
           
            persistent iteration
            persistent costlines
            persistent cstrlines
            if nargin == 1
               clear iteration;
               clear costlines;
               clear cstrlines;
               return;
            end
            % Initialize iteration count
            if isempty(iteration)
               iteration = 1; 
            end
            % Print counter, costs and constraints
            [cost, cstr] = obj.calculateCostsAndConstraints(z);
            costNames = fieldnames(cost);
            cstrNames = fieldnames(cstr);
            valNames = [costNames; cstrNames];
            args = cell(1,length(valNames));
            
            % Create the string formats to print out the costs
            titleStr = cell(1,length(valNames));
            formatStr = cell(1,length(valNames));

            for n = 1:length(valNames)
                titleStr{n} = [' ',valNames{n},' \t'];
                formatStr{n} = [' %',num2str(length(valNames{n})), '.6e \t'];
                if n <= length(costNames)
                    args{n} = sum(cost.(valNames{n}));
                else
                    args{n} = sum(cstr.(valNames{n}));
                end
            end
            
            titleStr = [' Iteration \t', cat(2,titleStr{:})];
            formatStr = [' %12d \t', cat(2,formatStr{:})];
            
            % Every 100 iterations, print the labels
            if rem(iteration, 100) == 1
               fprintf([titleStr, '\n']);
            end
            % Print the cost values
            fprintf([formatStr, '\n'], iteration, args{:});
            
            
%             % Make a figure
%             if iteration  == 1
%                 % Plot the figures
%                 costlines = figure();
%                 for n = 1:length(costNames)
%                    subplot(length(costNames),1, n);
%                    plot(iteration, sum(cost.(costNames{n})),'k');
%                    ylabel(costNames{n});
%                 end
%                 xlabel('Iterations');
%                 
%                 % Plot the constraints
%                 cstrlines = figure();
%                 for n = 1:length(cstrNames)
%                     subplot(length(cstrNames), 1, n)
%                     plot(iteration, sum(cstr.(cstrNames{n})), 'k');
%                     ylabel(cstrNames{n})
%                 end
%                 xlabel('Iterations');
%             else
%                 figure(costlines)
%                 for n = 1:length(costNames)
%                    subplot(length(costNames), 1, n)
%                    line = get(gca,'Children');
%                    set(line, 'XData', [get(line, 'XData'), iteration], 'YData', [get(line, 'YData'), sum(cost.(costNames{n}))]); 
%                 end
%                 figure(cstrlines)
%                 for n = 1:length(cstrNames)
%                    subplot(length(cstrNames), 1, n)
%                    line = get(gca,'Children');
%                    set(line,'XData', [get(line,'XData'), iteration], 'YData', [get(line,'YData'), sum(cstr.(cstrNames{n}))]);
%                 end
%             end
            % Increment iteration count
            iteration = iteration + 1;
            
        end
        function [costVals, cstrViol, cstrVal] = calculateCostsAndConstraints(obj,z)
            
            
            %% Cost function Values
            numCosts = length(obj.cost);
           
            costs = zeros(1,numCosts);
            costNames = cell(1,numCosts);
            % Calculate the values of the cost functions
            for n = 1:numCosts
                % Get the variable indices
                inds = obj.cost_xind_cell{n};
                % Get the variables
                args = cell(1,length(inds));
                for k = 1:length(inds)
                   args{k} = z(inds{k}); 
                end
                % Calculate the cost
                costs(n) = obj.cost{n}.eval_handle(args{:});
                % Store the costs and the labels
                costNames{n} = obj.cost{n}.name{1};
            end
            % Sort the costs by name
            [uniqueNames, ~, id] = unique(costNames);
            costVals = struct();
            % Sum the common costs
            for n = 1:length(uniqueNames)
               costVals.(uniqueNames{n}) = costs(id == n); 
            end
            
            %% Constraint Values
            numCstr = length(obj.nlcon);
            val = cell(1,numCstr);
            viols = zeros(1,numCstr);
            cstrName = cell(1,numCstr);
            % Calculate the values of the constraints
            for n = 1:length(obj.nlcon)
               % Get the variable indices
               inds = obj.nlcon_xind{n};
               % Get the variables
               args = cell(1,length(inds));
               for k = 1:length(inds)
                  args{k} = z(inds{k}); 
               end
               % Now get the value of the constraint
               val{n} = obj.nlcon{n}.eval_handle(args{:});
               % Calculate the maximum violation
               viol = max(obj.nlcon{n}.lb - val{n}, val{n} - obj.nlcon{n}.ub);
               % Eliminate the infinities
               viol(isinf(viol)) = 0;
               % Take the norm
               viols(n) = norm(viol);
               % Add to the associated labeled constraint
               cstrName{n} = obj.nlcon{n}.name{1};
            end
%             for n = 1:length(obj.bbcon)
%                 % Get the arguments to the bounding box constraints
%                 k = length(obj.nlcon) + n;
%                 inds = obj.bbcon_xind{n};
%                 args = z(inds);
%                 % Get the value of the constraints
%                 val{k} = obj.bbcon{n}.A * args;
%                 viol = max(obj.bbcon{n}.lb - val{k}, val{k} - obj.bbcon{n}.ub);
%                 viol(isinf(viol)) = 0;
%                 viols(n) = norm(viol);
%                 cstrName{k} = obj.bbcon{n}.name{1};
%             end
            % Sum the common constraints
            [uniqueNames, ~, id] = unique(cstrName);
            cstrVal = struct();
            cstrViol = struct();
            for n = 1:length(uniqueNames)
                cstrVal.(uniqueNames{n}) = cat(2,val{id == n});
                cstrViol.(uniqueNames{n}) = viols(id == n);
            end
        end
        function cstrVals = calculatContactConstraints(obj, z)
            %% calculateContactConstraints: Helper function for calculating contact constraints
            cstrVals = zeros(obj.numContacts * (2 + obj.numFriction), obj.N-1);
            for n = 1:obj.N-1
               cstrVals(:,n) =  obj.contact_constraint_fun(z([obj.x_idns(:,n+1); obj.lambda_inds(:,n)]));
            end
        end
    end
    methods (Access = protected)
        function [f, df] = relaxed_nc_cost(obj, x)
            % scaledSumCost: 
            f = obj.options.relax_cost * sum(x(:));
            df = obj.options.relax_cost * ones(1, numel(x));
        end
    end
end