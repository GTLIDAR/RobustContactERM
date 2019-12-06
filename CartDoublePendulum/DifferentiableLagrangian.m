classdef DifferentiableLagrangian 
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        numQ (1,1) double {mustBePositive} = 1;
        numU (1,1) double {mustBeNonnegative} = 1;
    end
    
    methods (Abstract)
       [w, dw, ddw] = configurationTransformation(obj, q);
       [B, dB] = controllerMatrix(obj, q);
    end
    
    properties (Abstract, Dependent)
       massFunctional;
       potentialFunctional;
       positionFunctional;
    end
    
    methods
        function obj = DifferentiableLagrangian(nQ, nU)

           obj.numQ = nQ;
           obj.numU = nU;
        end
        function [f, df] = dynamics(obj, ~, x, u)
            %% DYNAMICS: Evaluates the dynamics of the Lagrangian System
            %
            %   [f, df] = dynamics(plant, t, x, u) returns the dynamics of the
            %   Lagrangian System PLANT as a vector of derivatives of the
            %   state such satisfying:
            %           f = dx/dt
            %
            %   In general, we write the dynamics as a first-order form as:
            %           dx/dt = f(x,u)
            %
            %   Where x is the state vector and u is the control vector
            %
            %   DYNAMICS also returns the gradient of the dynamics, df,
            %   with respect to both the state and the controls:
            %   
            %       df = [df/dx; df/du]
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
            
            % Get the configuration and the velocities
            q = x(1:obj.numQ);
            dq = x(obj.numQ+1:end);
            
            % Get the physical parameters of the system (mass matrix, etc)
            %[M, dM] = obj.massMatrix(q);
            [C, dC, M, dM] = obj.coriolisMatrix(q,dq);  %Optionally return the mass matrix, to avoid calculating it twice.
            [N, dN] = obj.gravityMatrix(q);
            [B, dB] = obj.controllerMatrix(q);

            % Invert the mass matrix, as we'll use the inverse several
            % times.
            Minv = M\eye(obj.numQ);
                        
            % Calculate the dynamics f
            tau = B*u - C*dq - N;
            ddq = Minv*tau;
            f = [dq;ddq];
            if nargout == 2
                % Calculate the gradient wrt q
                dBu = times(dB, reshape(u, 1, numel(u), 1));
                dBu = squeeze(sum(dBu,2));
                dCv = times(dC(:,:,1:obj.numQ), reshape(dq, 1, obj.numQ, 1));
                dCv = squeeze(sum(dCv, 2));
                
                % Gradient of tau wrt q
                dtau = dBu - dCv - dN;                         
                
                % Gradient of the inverse mass matrix
                dMinv = zeros(size(dM));
                for n = 1:obj.numQ
                    dMinv(:,:,n) = -Minv * dM(:,:,n) * Minv;
                end
                % Multiply by tau
                dMinv_tau = times(dMinv, reshape(tau, 1, obj.numQ, 1));
                dMinv_tau = squeeze(sum(dMinv_tau, 2));
                
                % Gradient of the dynamics wrt q
                dfdq = [zeros(obj.numQ), dMinv_tau + Minv*dtau;];
                
                % Gradient wrt dq
                dCv_dq = times(dC(:,:,obj.numQ+1:end), reshape(dq, 1, obj.numQ, 1));
                dCv_dq = squeeze(sum(dCv_dq,2));
                dfddq = [eye(obj.numQ); M\(dCv_dq + C);];
                
                % Gradient wrt u
                dfdu =[zeros(obj.numQ); M\B];
                
                % Combine all the gradients into one
                df = [zeros(2*obj.numQ, 1), dfdq, dfddq, dfdu];
            end
        end
        function [H, C, B, dH, dC, dB] = manipulatorDynamics(obj, q, dq)
            %% MANIPULATORDYNAMICS: Returns mass matrix and forces for calculating manipulator dynamics
            %
            %   [H, C, B, dH, dC, dB] = manipulatorDynamics(PLANT, q, dq)
            %   calculates mass matrix H and force effects C for the PLANT
            %   to form the second-order equations of motion:
            %
            %       H*ddq + C = B*u
            %
            %   where u is the control input and B is the control selection
            %   matrix. MANIPULATORDYNAMICS also returns the gradients of
            %   the mass matrix dH, the force effects dC, and the control
            %   selection matrix dB.
            %
            %   MANIPULATORDYNAMICS is required to implement the Drake
            %   Class MANIPULATOR.
            %
            %   ARGUMENTS:
            %       PLANT:  a DIFFERENTIABLELAGRANGIAN Model
            %       q:      nQ x 1 vector of configuration variables
            %       dq:     nQ x 1 vector of configuration rates
            %
            %   RETURN VALUES:
            %       H:      nQ x nQ matrix, the mass matrix
            %       C:      nQ x 1 vector of coriolis and gravitational
            %               effects
            %       B:      nQ x nU matrix, the controller selection matrix
            %       dH:     nQ^2 x nQ matrix of partial derivatives of the
            %               mass matrix. dH=reshape(dH', [nQ, nQ, nQ])
            %               reorganizes the matrix such that dH(:,:,k) is
            %               the partial derivative of H wrt q(k)
            %       dC:     2*nQ x nQ matrix partial derivatives of force
            %               effects. The first nQ rows are the partial 
            %               derivatives wrt q, and the second nQ rows
            %               are the partial derivatives wrt dq
            %       dB:     nU*nQ x nQ matrix of partial derivatives of the
            %               controller selection matrix B.
            %               dB = reshape(dB', [nQ, nU, nQ]) reorganizes the
            %               array such that the dB(:,:,k) is the partial 
            %               derivative with respect to q(k). 
            
            % Get the dynamic properties
            [C,dC,H,dH] = obj.coriolisMatrix(obj,q,dq);
            [N, dN] = obj.gravityMatrix(q);
            [B, dB] = obj.controllerMatrix(q);
            
            % Gradient wrt q
            dC_q = times(dC(:,:,1:obj.numQ), reshape(dq,1,obj.numQ,1));
            dC_q = squeeze(sum(dC_q, 2));
            
            dC_q = dC_q + dN;
            % Gradient wrt dq
            dC_dq = C;
            
            % Gradient of C*q + N
            dC = [dC_q, dC_dq]';
            
            % Combine the coriolis and gravitational effects
            C = C*dq + N;
            
            % Reshape dM and dB
            dH = reshape(dH, [size(dH, 1), size(dH,2) * size(dH, 3)])';
            dB = reshape(dB, [size(dB, 1), size(dB, 2) * size(dB, 3)])';
        end
        function [E, T, V] = energy(obj, q, dq)
            %% ENERGY: Total Mechanical Energy of the Plant
            %
            %   [E, T, V] = energy(plant, q, dq) returns the total
            %   mechanical energy E of the PLANT at configuration Q and 
            %   rate DQ. ENERGY also optionally returns the kinetic energy
            %   T and the potential energy V.
            %
            %   ARGUMENTS:
            %       PLANT:  a DIFFERENTIABLELAGRANGIAN Model
            %       q:      nQ x 1 vector of configuration variables
            %       dq:     nQ x 1 vector of configuration rates
            %
            %   RETURN VALUES:
            %       E:      scalar double, the total energy. E = T + V
            %       T:      scalar double, the kinetic energy
            %       V:      scalar double, the potential energy
            
            
            % Calculate the Kinetic Energy
            M = obj.massMatrix(q);
            T = 0.5 * dq' * M * dq;
            % Calculate the Potential Energy
            v = obj.potentialFunctional;
            w = obj.configurationTranformation(q);
            V = v * w;
            % Total Energy is Kinetic + Potential            %   Arguments;
            E = T + V;
        end
        function [M, dM, ddM] = massMatrix(obj, q)
            %% MASSMATRIX: Matrix of Intertial Effects
            %
            %   [M, dM, ddM] = massMatrix(plant, q) returns the matrix of
            %   inertial effects M, its gradient, dM, and its Hessian ddM,
            %   for the PLANT model at configuration q. The Mass Matrix is 
            %   used in the calculation of the PLANT dynamics:
            %
            %       M*ddq + C*dq + N = B*u
            %
            %   The gradient dM and the Hessian ddM are used in the
            %   calculation of the gradients of the dynamics and in the
            %   calculation of the Corilois Matrix and its gradient.
            %
            %   ARGUMENTS:
            %       PLANT:  a DIFFERENTIABLELAGRANGIAN Model
            %       q:      nQ x 1 vector of configuration variables
            %
            %   RETURN VALUES:
            %       M:      nQ x nQ matrix of inertial effects
            %       dM:     nQ x nQ x nQ array of gradients for M.
            %               dM(:,:,k) is the derivative of M with respect to
            %               q(k)
            %       ddM:    nQ x nQ x nQ x nQ array of Hessians (second
            %               derivatives). ddM(:,:,i,j) is the mixed partial
            %               derivative if M with respect to q(i) and q(j).

            % Get the transformation and the derivative vectors.
            [w, dw, ddw] = obj.configurationTransformation(q);
            % The mass functional matrix
            Mbar = obj.massFunctional;
            % Duplicate to get the full mass matrix
            DMbar = duplication(obj.numQ) * Mbar;
            % Form the mass matrix
            M = reshape(DMbar * w, obj.numQ * ones(1,2));
            % Form the gradient
            dM = reshape(DMbar * dw, obj.numQ * ones(1,3));  
            % Form the hessian
            ddM = reshape(DMbar * ddw, obj.numQ * ones(1,4)); 
        end            
        function [C, dC, M, dM] = coriolisMatrix(obj, q, dq)
            %% CORIOLISMATRIX: Matrix of Coriolis and Centripetal Effects
            %
            %   [C,dC] = coriolisMatrix(plant, q, dq) returns the matrix of
            %   Coriolis and Centripetal Effects, C, and its gradient, dC,
            %   for the PLANT model at configuration q and configuration
            %   rate dq. The Coriolis Matrix is used in the calculation of
            %   the PLANT dynamics:
            %
            %       M*ddq + C*dq + N = B*u
            %
            %   ARGUMENTS:
            %       PLANT:  a DIFFERENTIABLELAGRANGIAN Model
            %       q:      nQ x 1 vector of configuration variables
            %       dq:     nQ x 1 vector of configuration rates
            %
            %   RETURN VALUES:
            %       C:      nQ x nQ double matrix of Coriolis and
            %           centripetal effects
            %       dC:     nQ x nQ x 2nQ array of gradients for the 
            %           Coriolis matrix C. The first nQ pages are the
            %           gradients with respect to q, the second nQ pages
            %           are the gradients with respect to dq
            
            % Get the mass matrix gradient and Hessian to calculate the Coriolis Matrix
            [M, dM, d2M] = obj.massMatrix(q);
            % Calculate the chirstoffel symbols
            G = 0.5 * (dM + permute(dM, [1,3,2]) - permute(dM, [3,2,1]));
            % Calculate the partials of the Christoffel symbols
            dG = 0.5 * (d2M + permute(d2M,[1,3,2,4]) - permute(d2M, [3,2,1,4]));
            % Calculate the Coriolis Matrix and its gradient using tensor multiplication
            C = times(G,reshape(dq,1,1,obj.numQ));
            C = squeeze(sum(C,3));
            dG = times(dG, reshape(dq,1,1,obj.numQ,1));
            dC = squeeze(sum(dG,3));
            % The gradient wrt dq is the Christoffel Symbols, G
            dC = cat(3, dC, G);             
            
        end
        function [N, dN] = gravityMatrix(obj, q)
            %% GRAVITYMATRIX: Matrix of Gravitational and Conservative Forces
            %
            %   [N, dN] = gravityMatrix(plant, q) returns the vector of
            %   gravitational and conservative Effects, N, and its gradient
            %   dN for the PLANT model at configuration q. The conservative
            %   force vector N is used in the calculation of the dynamics:
            %
            %       M*ddq + C*dq + N = B*u
            %
            %   ARGUMENTS:
            %       PLANT:  a DIFFERENTIABLELAGRANGIAN Model
            %       q:      nQ x 1 vector of configuration variables
            %
            %   RETURN VALUES:
            %       N:      nQ x 1 vector of conservative forces
            %       dN:     nQ x nQ matrix of gradients for the 
            %               conservative forces with respect to q.
            
            % Get the potential energy functional 
            v = obj.potentialFunctional;   
            [~, dw, ddw] = obj.configurationTransformation(q);
            % Calculate the conservative forces and their gradient
            N = (v * dw)';
            dN = reshape((v * ddw)', obj.numQ, obj.numQ);
        end
        function [J, dJ] = jacobian(obj, q)
            %% Jacobian: Returns the Jacobian mapping configuration rates to Cartesian Velocities
            %
            %   [J, dJ] = jacobian(plant, q) returns the Jacobian matrix J 
            %       and its gradient dJ for the PLANT model evaluated at 
            %       the configuration q.
            %
            %   ARGUMENTS:
            %       PLANT:  a DIFFERENTIABLELAGRANGIAN Model
            %       q:      nQ x 1 vector of configuration variables
            %
            %   RETURN VALUES:
            %       J:      D x nQ Jacobian matrix
            %       dJ:     D x nQ x nQ array of gradients for the Jacobian 
            %               matrix. dJ(:,:,k) is the partial derivative of
            %               J with respect to q(k).
            
            % Get the configuration transformation derivatives and the
            % position functional
            [~, dw, ddw] = obj.configurationTransformation(q);
            P = obj.positionFunctional;
            % The Jacobian and its gradient
            J = P * dw;
            dJ = reshape(P * ddw, size(P,1), obj.numQ, obj.numQ);
        end
        function x = kinematics(obj, q)
            %% Kinematics: Calculates the Cartesian coordinates of the endpoint of the system
            %
            %   x = kinematics(plant, q) returns the Cartesian coordinates x
            %   of the endpoint of PLANT at configuration q.
            %
            %   ARGUMENTS:
            %       PLANT:  a DIFFERENTIABLELAGRANGIAN Model
            %       q:      nQ x 1 vector of configuration variables
            %
            %   RETURN VALUES:
            %       x:      D x 1 vector of Cartesian coordinates

            w = obj.configurationTransformation(q);
            P = obj.positionFunctional;            
            x = P*w;
        end
    end
end

