classdef DifferentiableContactDynamics 
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        terrain;
        timestep = 0.01;
        contactSolver = PathLCPSolver();  
    end
    
    methods
        
        function obj = DifferentiableContactDynamics(terrainObj)
           if nargin == 1
               obj.terrain = terrainObj;
           else
              obj.terrain = FlatTerrain(); 
           end
        end
    end
    methods (Sealed = true)
        function [t,x] = simulate(obj, Tf, x0)
            
            % Pre-initialize arrays
            nX = numel(x0);
            t = 0:obj.timestep:Tf;
            nT = numel(t);
            x = zeros(nX, nT);
            x(:,1) = x0;
            u = zeros(obj.numU, 1);
            % Run the simulation loop
            for n = 1:nT-1
               dx = obj.dynamics(t(n), x(:,n), u);
               x(:,n+1) = x(:,n) + dx * obj.timestep;
            end
        end
        function [f, df] = dynamics(obj, t, x, u)
            %% DYNAMICS: Evaluates the dynamics of the Lagrangian System
            %
            %   [f, df] = dynamics(plant, t, x, u) returns the dynamics of the
            %   Lagrangian System PLANT as a vector of derivatives of the
            %       xA: the point on the manipulator nearest the terrain
            %   state such satisfying:
            %           f = dx/dt
            %
            %   In general, we write the dynamics as a first-order form as:
            %           dx/dt = f(t,x,u)
            %
            %   Where x is the state vector and u is the control vector
            %
            %   DYNAMICS also returns the gradient of the dynamics, df,
            %   with respect to both the state and the controls:
            %   
            %       df = [df/dt; df/dx; df/du]
            %      
            %   ARGUMENTS:
            %       PLANT:  a DIFFERENTIABLELAGRANGIAN Model
            %       t:      scalar time variable (unused)
            %       x:      2*nQ x 1 vector of state variables
            %       u:      2*nQ x 1 vector of control variables
            %
            %   RETURN VALUES:
            %       f:      2*nQx1 vector of first derivatives of the state
            %       df:     2*nQ x (2*nQ + nU + 1) vector, the gradient of the
            %               dynamics
            %
            %   The returned dynamics are evaluated semi-implicitly. That
            %   is, they are evaluated at the next time step, assuming
            %   constant velocity.
            
            % Get the configuration and the velocities
            q = x(1:obj.numQ);
            dq = x(obj.numQ+1:end);
            % Solve the contact problem
            [fc, dfc, ~] = obj.contactForce(q, dq, u);
            
            fc = fc./obj.timestep;
            dfc = dfc./obj.timestep;
            
            qhat = q + obj.timestep * dq;
            % Get the physical parameters of the system (mass matrix, etc)
            [M, dM] = obj.massMatrix(qhat);
            [C, dC] = obj.coriolisMatrix(qhat,dq);  %Optionally return the mass matrix, to avoid calculating it twice.
            [N, dN] = obj.gravityMatrix(qhat);
            [B, dB] = obj.controllerMatrix(qhat);
            [Jn, Jt, dJn, dJt] = obj.contactJacobian(qhat);
            % Transpose the Jacobians
            Jc = [Jn; Jt]';
            dJc = permute(cat(1,dJn, dJt), [2,1,3]);
            
            % Invert the mass matrix, as we'll use the inverse several
            % times.
            Minv = M\eye(obj.numQ);
                        
            % Calculate the dynamics f
            tau = B*u - C*dq - N;
            ddq = Minv*(tau + (Jc * fc)./obj.timestep);
            
            % Time-stepping dynamics (semi-implicit, evaluated at the next time
            % step)
            f = [dq + obj.timestep*ddq; ddq];
            
            if nargout == 2
                % Split dfc into parts corresponding to q, dq, and u
                dfc_q = dfc(:,1:obj.numQ);
                dfc_dq = dfc(:,obj.numQ+1:2*obj.numQ);
                dfc_u = dfc(:,2*obj.numQ+1:end);
                
                % Calculate the gradient wrt q
                dBu = times(dB, reshape(u, [1, numel(u), 1]));
                dBu = squeeze(sum(dBu,2));
                dC_dq = times(dC, reshape(dq, [1, obj.numQ, 1]));
                dC_dq = squeeze(sum(dC_dq, 2));
                
                dJc_f = times(dJc, reshape(fc, [1, numel(fc), 1]));
                dJc_f = squeeze(sum(dJc_f, 2));
                
                % Gradient of tau wrt q
                dtau_q = dBu - dC_dq(:,1:obj.numQ) - dN;% + dJc_f + Jc * dfc_q;                         
                
                % Gradient of tau wrt dq
                dtau_dq = obj.timestep * dtau_q - (dC_dq(:,obj.numQ+1:end) + C);
                
                % Gradient of the inverse mass matrix
                dMinv = zeros(size(dM));
                for n = 1:obj.numQ
                    dMinv(:,:,n) = -Minv * dM(:,:,n) * Minv;
                end
                % Multiply by tau
                dMinv_tau = times(dMinv, reshape(tau, [1, obj.numQ, 1]));
                dMinv_tau = squeeze(sum(dMinv_tau, 2));            
                % Multiply by fc
                dMinv_fc = times(dMinv, reshape(Jc*fc, [1, obj.numQ, 1]));
                dMinv_fc = squeeze(sum(dMinv_fc,2));
                        
                % Calculate terms common to df2_q, df2_dq
                df2_common = dMinv_tau + (dMinv_fc + Minv*dJc_f)/obj.timestep;
                
                % Calculate the partial derivatives
                df2_q = df2_common + Minv*dtau_q + Minv*Jc*dfc_q / obj.timestep;
                
                df2_dq = obj.timestep * df2_common + Minv*(dtau_dq + Jc*dfc_dq/obj.timestep);
                
                df2_u = Minv * B + Minv * (Jc * dfc_u)./obj.timestep;
                
                df_q = [obj.timestep * df2_q; df2_q];
                df_dq = [eye(obj.numQ) + obj.timestep * df2_dq; df2_dq];
                df_u = [obj.timestep * df2_u; df2_u];
                                
                % Combine all the gradients into one
                df = [zeros(2*obj.numQ, 1), df_q, df_dq, df_u];
            end
        end 
        function [f,df,r] = contactForce(obj,q, dq,u)

            % Get the parameters of the LCP problem
            [P, z, dP, dz, numF] = obj.getLCP(q, dq, u);
            % Solve the LCP
            [f,r] = obj.contactSolver.solve(P,z);
            % The gradient of the force
            if any(f(1:numF))
                df = obj.contactSolver.gradient(f, P, z, dP, dz);
                % Get only the normal and tangential components
                df = df(1:numF,:);
            else
                df = zeros(numF, 2*obj.numQ + obj.numU);
            end
            % Get the normal and tangential forces
            f = f(1:numF,:);
        end

        function [P, z, dP, dz, numF] = getLCP(obj, q, dq, u)
            % Simulate the configuration using constant velocity
            qhat = q + obj.timestep*dq;
            % Get the system properties
            %[M, C, B, dM, dC, dB]  = obj.manipulatorDynamics(qhat, dq, u);
            [M, dM] = obj.massMatrix(qhat);
            [C, dC] = obj.coriolisMatrix(qhat,dq);
            [N, dN] = obj.gravityMatrix(qhat);
            [B, dB] = obj.controllerMatrix(qhat);
            
            % Invert the mass matrix, as we'll be using the inverse often
            iM = M\eye(obj.numQ);
            [Jn, Jt, dJn, dJt, alphas] = obj.contactJacobian(qhat);

            % Collect the contact Jacobian into one term:
            Jc = [Jn; Jt];
            dJc = cat(1,dJn, dJt);
            
            % Pre-initialize the arrays
            numT = size(Jt,1);
            numN = size(Jn,1);
            z = zeros(numT + 2*numN, 1);
            % Calculate the force effects, tau, and the no-contact velocity
            % vel:
            tau = B*u - C*dq - N;
            vel = dq + obj.timestep * iM * tau;
           
            % The LCP offset vector
            z(1:numN) = Jn*vel + (Jn*q - alphas)/obj.timestep;
            z(numN+1:numN+numT) = Jt*vel;

            % Calculate the problem matrix 
            P = zeros(numT+2*numN, numT+2*numN);
            w = zeros(numN,numT);
            
            for n = 1:numN
                w(n,numN*(n-1)+1:numN*n) = 1;
            end
                        
            P(1:numN + numT, 1:numN + numT) = Jc * iM * Jc';
            P(numN+1:numN+numT, numN+numT+1:end) = w';
            P(numN + numT+1:end,1:numN) = eye(numN) * obj.terrain.friction_coeff;
            P(numN+numT+1:end,numN+1:numN+numT) = -w;
            
            %% Gradient of the LCP Problem
            % Gradient of the problem matrix, dP
            % We'll also need the partial derivative of the inverse
            % mass matrix, so we'll calculate that here too.
            diM = zeros(obj.numQ*[1, 1, 1]);
            dP = zeros(numT+2*numN, numT+2*numN, 2*obj.numQ + obj.numU);
            for n = 1:obj.numQ
                diM(:,:,n) = -iM * dM(:,:,n) * iM;
                dP(1:numN + numT, 1:numN + numT,n) = dJc(:,:,n) * iM * Jc' + Jc * diM(:,:,n) * Jc' + Jc * iM * dJc(:,:,n)';
            end
            % Gradient of P wrt dq
            dP(:,:,obj.numQ+1:2*obj.numQ) = obj.timestep * dP(:,:,1:obj.numQ);
            % Calculate the gradient of the offset vector, z:
            % First we do some multidimensional array multiplication:
            dC = times(dC,reshape(dq, [1, obj.numQ, 1]));
            dC = squeeze(sum(dC,2));
            
            dBu = times(dB, reshape(u, [1, obj.numU, 1]));
            dBu = squeeze(sum(dBu, 2));
            
            % Gradients of TAU wrt q, dq
            dtau_q = dBu - dC(:,1:obj.numQ) - dN;
            dtau_dq = -(dC(:,obj.numQ+1:end) + C) + obj.timestep*dtau_q;
            
            % Tensor Multiplications
            diM_tau = times(diM, reshape(tau, [1, obj.numQ, 1]));
            diM_tau = squeeze(sum(diM_tau, 2));

            dJc_v = times(dJc, reshape(vel, [1, numel(vel), 1]));
            dJc_v = squeeze(sum(dJc_v, 2));
            
            % Common terms to dz_q, dz_dq
            dz_common = dJc_v + obj.timestep * Jc * diM_tau;
            
            % Gradient wrt q, dq
            dz_q = dz_common + obj.timestep* Jc * iM * dtau_q;
            dz_q(1:numN,:) = dz_q(1:numN,:) + Jn ./ obj.timestep;
            dz_dq = obj.timestep * dz_common + Jc * (eye(obj.numQ) + obj.timestep * iM * dtau_dq);
            dz_u = obj.timestep * Jc * iM * B;
            
            
            % Collect all the gradients in one
            dz = zeros(numT + 2*numN, 2*obj.numQ + obj.numU);
            dz(1:numT+numN, :) = [dz_q, dz_dq, dz_u];
            
            numF = numT + numN;
        end
        
        function [Jn, Jt, dJn, dJt, alphas] = contactJacobian(obj, q)
            
           % Find the nearest point on the terrain 
           xA = obj.kinematics(q);
           xB = obj.terrain.nearest(xA);
           
           % Calculate the gradient of the terrain at that point
           [N, T] = obj.terrain.basis(xB);
           T = [T, -T];
           % Get the Jacobian
           [J, dJ] = obj.jacobian(q);
           
           % Calculate the contact Jacobian and its derivatives
           Jn = N'*J;
           Jt = T'*J;
           
           dJn = zeros([size(Jn), obj.numQ]);
           dJt = zeros([size(Jt), obj.numQ]);
           
           for n = 1:size(dJ, 3)
              dJn(:,:,n) = N'*dJ(:,:,n);
              dJt(:,:,n) = T'*dJ(:,:,n);
           end
           
           % Calculate the alpha values
           phi = N' * (xA - xB);
           alphas = Jn * q - phi;
        end
                
        function [phi, Nw, Dw, xA, xB, iA, iB, mu, N, D, dN, dD] = contactConstraints(obj, q, ~,~) 
            %   Return Values
            %       phi: the signed distance between the manipulator and
            %       the ground
            %       Nw:  the contact normal expressed in world coordinates
            %       Dw: the friction basis expressed in world coordinates
            %       xA: the point on the manipulator nearest the terrain
            %       xB: the point on the terrain nearest the manipulator
            %       iA: the body index for the contact point xA (1)
            %       iB: the body index for the contact point xB (0)
            %       mu: the coefficient of friction
            %       N: the contact normal expressed in generalized
            %       coordinates
            %       D: the friction basis expressed in generalized
            %       coordinates
            %       dN: the derivative of N, dN/dq
            %       dD: the derivative of D, dD/dq
            
            % The body indices for the contactDrivenCart are fixed
            % quantities
            iA = 1; % Body A is the manipulator (the only body in the problem)
            iB = 0; % Body B is the terrain
            
            % Friction is a property of the terrian
            mu = obj.terrain.friction_coeff;
            
            % Get the nearest point from the terrain
            xA = obj.kinematics(q);
            xB = obj.terrain.nearest(xA);
            
            % Calculate the local coordinates of the terrain
            [Nw, Tw] = obj.terrain.basis(xB);
            Dw = {Tw};
            % Calculate thenumT + 2*numN, 1 signed distance function
            phi = Nw' * (xA - xB);
            phi = phi';
            % Calculate N and D from the jacobian
            [J, dJ] = obj.jacobian(q);
            N = J' * Nw;
            D = cell(2,1);
            D{1} = J' * Tw;
            D{2} = -D{1};
            % Calculate dN and dD by contraction with dJacobian
            dN = zeros(numel(N), numel(q));
            dD = cell(2);
            dD{1} = zeros(numel(D{1}),numel(q));
            dD{2} = dD{1};
            for n = 1:numel(q)
                tempN = dJ(:,:,n)'*Nw;
                dN(:,n) = tempN(:);
                tempD = dJ(:,:,n)'*Tw;
                dD{1}(:,n) = tempD(:);
                dD{2}(:,n) = -tempD(:);
            end
            % Format for output
            N = N';
            D{1} = D{1}';
            D{2} = D{2}';
        end
    end
    
end

