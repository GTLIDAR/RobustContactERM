classdef SemiImplicitTrajectoryOptimization < DirtranTrajectoryOptimization
   %% SemiImplicitTrajectoryOptimization: A wrapper class for DirtranTrajectoryOptimization
   %
   %    SemiImplicitTrajectoryOptimization is a wrapper class for
   %    DirtranTrajectoryOptimization that adds the use of semi-implicit
   %    Euler integration of system dynamics. It is assumed that the
   %    state can be partitioned as x = [x1, x2], where x1 and x2 have
   %    equal lengths, such that the dynamics:
   %        dx1 = x2
   %        dx2 = f(t, x, u);
   %    Such dynamics are common in Langrangian Systems. In this case,
   %    SemiImplicitTrajectoryOptimization implements a semi-implicit Euler
   %    integration scheme where x2 is updated using forward euler
   %    integration and x1 is updated using backward euler integration:
   %        x2(k+1) = x2(k) + h * f
   %        x1(k+1) = x1(k) + h *x2(k+1)
   %    
   %    The semi-implicit method is also consistent with the time-stepping
   %    method for rigid body systems with rigid contact.
   
   %    Luke Drnach
   %    January 23, 2020
    
    
    properties (Constant)
        % ADDITIONAL INTEGRATION METHODS
        % Semi-Implicit does a forward pass on dq and a backward pass on q,
        % consistent with the semi-implicit formulation for the contact
        % force
        SEMI_IMPLICIT = 4;
    end
    
    methods
        function obj = SemiImplicitTrajectoryOptimization(plant,N,duration,options)
            %% SemiImplicitTrajectoryOptimization: Class constructor method
            %
            %   SemiImplicitTrajectoryOptimization is a wrapper for
            %   DirtranTrajectoryOptimization, except that the default
            %   integrator is now the Semi-Implicit integrator. 
            %
            %   Syntax:
            %       obj = SemiImplicitTrajectoryOptimization(plant, N,
            %       duration, options)
            %
            %   Arguments:
            %       plant: 
            
            if nargin < 4
                options = struct();
            end
            if ~isfield(options,'integration_method')
                options.integration_method = SemiImplicitTrajectoryOptimization.SEMI_IMPLICIT;
            end
            obj = obj@DirtranTrajectoryOptimization(plant,N,duration,options);
        end
        
        function obj = addDynamicConstraints(obj)
            %% ADDDYNAMICCONSTRAINTS: Adds the dynamics as constraints to the trajectory optimization problem
            %
            %   
            
            
            % Check if we're using the Semi-Implicit Method. Otherwise,
            % default to using DirtranTrajectoryOptimization
            if(obj.options.integration_method == SemiImplicitTrajectoryOptimization.SEMI_IMPLICIT)
                % Setup
                nX = obj.plant.getNumStates();
                nU = obj.plant.getNumInputs();
                N = obj.N;
                % Initialize cell arrays
                constraints = cell(N-1,1);
                dyn_inds = cell(N-1,1);
                % Create the function handle for the constraints
                n_vars = 2*nX + nU + 1;
                cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.forward_constraint_fun);
                % Add the constraints at every point
                for i=1:obj.N-1
                    dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
                    constraints{i} = cnstr;
                    obj = obj.addConstraint(constraints{i}, dyn_inds{i});
                end
            else
                % Use the methods in DirtranTrajectoryOptimization (the
                % parent)
                obj = addDynamicConstraints@DirtranTrajectoryOptimization(obj);
            end
        end
        
        function obj = addRunningCost(obj,running_cost_function)
            % Adds an integrated cost to all time steps, which is
            % numerical implementation specific (thus abstract)
            % this cost is assumed to be time-invariant
            % @param running_cost_function a function handle
            %  of the form running_cost_function(dt,x,u)
            
            if (obj.options.integration_method == SemiImplicitTrajectoryOptimization.SEMI_IMPLICIT)
                nX = obj.plant.getNumStates();
                nU = obj.plant.getNumInputs();
                running_cost = FunctionHandleObjective(1+nX+nU, running_cost_function);
                for i = 1:obj.N-1
                    inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.u_inds(:,i)};
                    obj = obj.addCost(running_cost,inds_i);
                end
            else
                obj = addRunningCost@DirtranTrajectoryOptimization(obj, running_cost_function);
            end
        end
        
    end
    methods
        function [f, df] = semiimplicit_constraint_fun(obj, h, x0, x1, u0)
            %% SEMI-IMPLICIT CONSTRAINT FUNCTION: calculate the defects from semi-implicit integration
            %
            %   returns:
            %       f: vector of defects
            %       df: gradients of the defects
            
            % Get the 'position' and 'velocity'
            nX = obj.plant.getNumStates();
            nQ = nX/2;
            q0 = x0(1:nQ);
            v0 = x0(nQ+1:end);
            q1 = x1(1:nQ);
            v1 = x1(nQ+1:end);
            
            %Calculate the dynamics
            [xdot, dxdot] = obj.plant.dynamics(0, x0, u0);
            vdot = xdot(nQ+1:end,:);
            dvdot = dxdot(nQ+1:end,:);
            
            % Forward integrate the velocity (the defect)
            fv = v1 - v0 - h * vdot;
            % Backward integrate the position (the defect)
            fx = q1 - q0 - h * v1;
            %f(1:nQ) = x1(1:nQ) - x0(1:nQ) - h * (x0(nQ+1:end) + h * vdot);
            % Calculate the gradients of the defects
            dfv = [-vdot, -h*dvdot(:,2:nQ+1), -(eye(nQ) + h*dvdot(:,nQ+2:2*nQ+1)),...
                zeros(nQ), eye(nQ), -h * dvdot(:,2*nQ+2 : end)];
            dfx = [-v1, -eye(nQ), zeros(nQ), eye(nQ), -h*eye(nQ), zeros(nQ, numel(u0))];
            %dfx = [-x0(nQ+1:end) - 2 * h * vdot, -(eye(nQ) + h^2 * dvdot(:,2:nQ+1)),...
                %-h * (eye(nQ) + h*dvdot(:,nQ+2:2*nQ+1)), eye(nQ), zeros(nQ), -h^2*dvdot(2*nQ+2:end)];
            % Combine
            f = [fx;fv];
            df = [dfx; dfv];
        end
    end
end