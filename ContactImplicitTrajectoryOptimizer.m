classdef ContactImplicitTrajectoryOptimizer < DirectTrajectoryOptimization
   
    properties
        numContacts;    % Number of contacts
        numFriction;    % Number of friction basis vectors per contact
        lambda_inds;    % Indices for finding the contact variables inside the decision variables, ordered [normal; tangential; velocity_slacks];
        force_inds;     % Indices for finding the force variables from inside the contact variables, ordered [normal; tangential];
    end
    properties (Constant)
        % INTEGRATION METHODS
        FORWARD_EULER = 1;
        BACKWARD_EULER = 2;
        MIDPOINT = 3;
        SEMI_IMPLICIT = 4;
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
               options.integration_method = ContactImplicitTrajectoryOptimizer.MIDPOINT; 
            end
            if ~isfield(options,'nlcc_mode')
                options.nlcc_mode = 1;
            end
            if ~isfield(options,'compl_slack')
                options.compl_slack = 0;
            end
            
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
           % Get the forces from the lambda variables
           obj.force_inds = 1:nContactForces;
           skip = 2 + obj.numFriction;
           vel_slacks = skip:skip:nContactForces;
           obj.force_inds(vel_slacks) = [];
           
           % Add decision variables for the forces
           obj = obj.addDecisionVariable(N * nContactForces);
        end
        function [xtraj, utraj, ftraj,z,F, info] = solveTraj(obj, t_init, traj_init)
            % Solve the problem using
            [xtraj, utraj, z, F, info] = solveTraj@DirectTrajectoryOptimization(obj, t_init, traj_init);
            
            % Pull the contact forces from the decision variable list
            if obj.numContacts > 0
                t = [0; cumsum(z(obj.h_inds))];
                ftraj = PPTrajectory(foh(t, reshape(z(obj.lambda_inds),[],obj.N)));
            else
                ftraj = [];
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
                    n_vars = 2*nX + nU + 1 + nL;
                    cstr = FunctionHandleConstraint(zeros(nX, 1), zeros(nX, 1), n_vars, @obj.forward_constraint_fun);
                case ContactImplicitTrajectoryOptimizer.BACKWARD_EULER
                    n_vars = 2*nX + nU + 1 + nL;
                    cstr = FunctionHandleConstraint(zeros(nX, 1), zeros(nX, 1), n_vars, @obj.backward_constraint_fun);
                case ContactImplicitTrajectoryOptimizer.MIDPOINT
                    n_vars = 2*nX + 2*nU + 1 + nL;
                    cstr = FunctionHandleConstraint(zeros(nX, 1), zeros(nX, 1), n_vars, @obj.midpoint_constraint_fun);
                case ContactImplicitTrajectoryOptimizer.SEMI_IMPLICIT
                    n_vars = 2*nX + nU + 1 + nL;
                    cstr = FunctionHandleConstraint(zeros(nX, 1), zeros(nX, 1), n_vars, @obj.semiimplicit_constraint_fun);
                otherwise
                    error('Unknown Integration Method');
            end
            % Get the indices necessary for the constraints (indices to
            % pull the states, controls, and forces from the decision
            % variable list).
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
            % Add in the complementarity constriants
            obj = obj.addContactConstraints();
        end   
        function obj = addContactConstraints(obj)
            nX = obj.plant.getNumStates();
            for i = 1:obj.N - 1
                % Add in the nonlinear complementarity constraints
                ncp_cstr = NonlinearComplementarityConstraint_original(@obj.contact_constraint_fun, nX, obj.numContacts * (2 + obj.numFriction), obj.options.nlcc_mode, obj.options.compl_slack);
                %ncp_cstr = NonlinearComplementarityConstraint(@contat_constraint_fun, nX, obj.numContacts*(2+obj.numFriction), obj.options.nlcc_mode);
                obj = obj.addConstraint(ncp_cstr, [obj.x_inds(:,i+1); obj.lambda_inds(:,i)]);
            end
        end
        
        function obj = addRunningCost(obj, running_cost_function)
           %% addRunningCost: add the running cost function as the objective
           %
           %    
            nX = obj.plant.getNumStates();
            nU = obj.plant.getNumInputs();
            
            for i=1:obj.N-1
                switch obj.options.integration_method
                    case ContactImplicitTrajectoryOptimizer.FORWARD_EULER
                        running_cost = FunctionHandleObjective(1+nX+nU, running_cost_function);
                        inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.u_inds(:,i)};
                    case ContactImplicitTrajectoryOptimizer.BACKWARD_EULER
                        running_cost = FunctionHandleObjective(1+nX+nU, running_cost_function);
                        inds_i = {obj.h_inds(i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
                    case ContactImplicitTrajectoryOptimizer.MIDPOINT
                        running_cost = FunctionHandleObjective(1+2*nX+2*nU,...
                            @(h,x0,x1,u0,u1) obj.midpoint_running_fun(running_cost_function,h,x0,x1,u0,u1));
                        inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1)};
                    case ContactImplicitTrajectoryOptimizer.SEMI_IMPLICIT
                        running_cost = FunctionHandleObjective(1+nX+nU, running_cost_function);
                        inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.u_inds(:,i)};
                    otherwise
                        error('Unknown integration method');
                end
                
                obj = obj.addCost(running_cost,inds_i);
            end
        end
    end
    methods 
        function [f, df] = forward_constraint_fun(obj, h, x0, x1, u, lambda)
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
           % Combine the defects and the derivatives
           f = [fq; fv];
           df = [dfq; dfv];
        end
        function [f, df] = backward_constraint_fun(obj, h, x0, x1, u, lambda)
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
            % Calculate the derivative
            dHv = squeeze(sum(dH .* (v1 - v0)', 2));
            dBu = squeeze(sum(dB .* u', 2));
            dJl = squeeze(sum(dJ .* lambda', 2));
            dfv = [-(B*u - C + J'*lambda), zeros(nV, nQ), -H, dHv - h*(dBu - dC(:,1:nQ) + dJl), H + h*dC(:,nQ+1:nQ+nV), -h*B, -h*J'];
            % Combine the defects and the derivatives
            f = [fq; fv];
            df = [dfq; dfv];
        end
        function [f, df] = midpoint_constraint_fun(obj, h, x0, x1, u0, u1, lambda)
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
            % Combine the defects and the derivatives
            f = [fq; fv];
            df = [dfq; dfv];
        end
        function [f, df] = semiimplicit_constraint_fun(obj, h, x0, x1, u, lambda)
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
            D = cat(3,D{:});
            dD = cat(3,dD{:});
            
            % Reshape the gradients
            dN = reshape(dN,[obj.numContacts, numel(q), numel(q)]);
            dD = reshape(dD,[obj.numContacts, numel(q), numel(q), obj.numFriction]);
            dD = permute(dD, [4,2,3,1]);
            
            for k = 1:obj.numContacts
                % Jacobian for the kth contact
                J{k} = [N(k,:); squeeze(D(k,:,:))'];
                % Jacobian derivative for the kth contact
                dJ{k} = cat(1,dN(k,:,:),dD(:,:,:,k));
            end
            J = cat(1, J{:});
            dJ = cat(1,dJ{:});
        end
        function [g, dg] = midpoint_running_fun(obj, cost_fun, h, x0, x1, u0, u1)
            
            nX = obj.plant.getNumStates();
            nU = obj.plant.getNumInputs();
            [g, dh] = cost_fun(h, 0.5*(x0 + x1), 0.5*(u0 + u1));
            
            dg = [dh(:,1), 0.5*dh(:,2:1 + nX), 0.5*dh(:, 2:1 + nX), 0.5 * dh(:, nX+2 : nX+nU+1), 0.5 * dh(:, nX+2: nX + nU + 1)];
        end
        function [f, df] = contact_constraint_fun(obj, y)
            %% CONTACT_NCP_CONSTRAINT_FUN: Function evaluating the contact nonlinear complementarity function
            %
            % CONTACT_NCP_CONSTRIANT_FUN implements the function whose values
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
            [phi, ~, ~, ~, ~, ~, ~, mu, N, D, ~, dD] = obj.plant.contactConstraints(q, false, obj.options.active_collision_options);
            % Reshape the tangential force basis for easier use
            D = cat(3,D{:});
            D = permute(D,[3,2,1]);
            dD = cat(3,dD{:});
            dD = reshape(dD, [obj.numContacts, nQ, nQ, size(dD, 3)]);
            dD = permute(dD, [4, 2, 3, 1]);
            % Initialize the complementarity function
            f = zeros(nZ, 1);
            df = zeros(nZ, nQ + nV + obj.numContacts*(2 + obj.numFriction));
            
            % Distance | Normal complementarity
            f(N_idx,:) = phi;       % Distance function
            df(N_idx,1:nQ) = N;     % Normal vector
                
            for i = 1:obj.numContacts
                % Slack | Cone complementarity
                f(G_idx(i),:) = lambda_N(i)*mu - sum(lambda_T(obj.numFriction*(i-1)+1 : obj.numFriction*i));
                % Derivative wrt normal force
                df(G_idx(i),nQ + nV + N_idx(i)) = mu;
                % Derivative wrt tangential force
                df(G_idx(i),nQ + nV + T_idx(obj.numFriction*(i-1)+1):nQ + nV + T_idx(obj.numFriction*i)) = -1;
                % Velocity | Tangential complementarity 
                f(T_idx(obj.numFriction*(i-1)+1):T_idx(obj.numFriction*i),:) = gamma(i) + D(:,:,i)*v;
                % Derivative wrt configuration
                df(T_idx(obj.numFriction*(i-1)+1):T_idx(obj.numFriction*i), 1:nQ) = squeeze(sum(dD(:,:,:,i).*v',2));
                % Derivative wrt velocity
                df(T_idx(obj.numFriction*(i-1)+1):T_idx(obj.numFriction*i), nQ+1:nQ+nV) = D(:,:,i);
                % Derivative wrt gamma
                df(T_idx(obj.numFriction*(i-1)+1):T_idx(obj.numFriction*i), nQ+nV+G_idx(i)) = 1;
            end
        end
    end
end