classdef ContactDynamics < LagrangianModel
    
    properties (Abstract)
       timestep; 
       friction;
    end
    properties
       contactSolver = "Path"; 
       sigma = 10;
    end
    properties (Hidden, Access=private)
       solver; 
       func;
       opts = optimoptions('lsqnonlin','MaxFunctionEvaluations',2000,'StepTolerance',1e-8);
       lastSoln = [];
    end
    methods (Abstract)
       [n,a] = contactNormal(self,q);
       J = contactJacobian(self,q);
    end
    methods
        function [time,X,F, r] = simulate(self,x0,T)
            self.checkSolver();
            % Pre-initialize arrays
            nX = numel(x0);
            time = 0:self.timestep:T;
            nT = numel(time);
            X = zeros(nX,nT);
            X(:,1) = x0;
            % Determine the size of the contact force
            [Jn, Jt] = self.contactJacobian(x0(1:nX/2));
            nF = size(Jt,2) + size(Jn,2);
            F = zeros(nF,nT);
            r = zeros(1,nT);
            % Run the main simulation loop
            for t = 1:nT-1
               % Calculate contact forces
               [F(:,t), r(t)] = self.contactForce(X(:,t));
               % Calculate the updated state vector
               X(:,t+1) = self.implicitDynamics(X(:,t),F(:,t));
            end
        end
        function x = implicitDynamics(self, x, f)
            % Get the configuration and generalized velocity
            L = length(x);
            q = x(1:L/2);
            qdot = x(L/2+1:end);
            % Estimate the next position, assuming constant velocity
            qhat = q + self.timestep*qdot;
            % Get the system properties
            M = self.inertiaMatrix(qhat);
            C = self.coriolisMatrix(qhat,qdot);
            N = self.gravityMatrix(qhat);
            [Jn, Jt] = self.contactJacobian(qhat);
            J = [Jn, Jt];
            % Update the generalized velocity
            %tau = self.timestep*(C*qdot + N) - M*qdot;
            %qdot = M\(self.timestep*J*f - tau);
            tau = self.timestep*(C*qdot +N);
            qdot = qdot + M\(J*f - tau);
            % Update the configuration
            q = q + self.timestep*qdot;
            % Return the state vector
            x = [q(:); qdot(:)];
        end
        function [f,r] = contactForce(self,x)
            % Use an LCP (a type of QP) to solve for the contact forces
            
            % Separate the configuration and velocity
            L = length(x);
            q = x(1:L/2);
            qdot = x(L/2+1:end);
            % Simulate the configuration using constant velocity
            qhat = q + self.timestep*qdot;
            % Get the system properties
            M = self.inertiaMatrix(qhat);
            C = self.coriolisMatrix(qhat,qdot);
            N = self.gravityMatrix(qhat);
            [Jn, Jt] = self.contactJacobian(qhat);
            [normals,alphas] = self.contactNormal(qhat);
            % Pre-initialize the arrays
            numT = size(Jt,2);
            numN = size(Jn,2);
            z = zeros(numT + 2*numN, 1);
            % Calculate the offset vector
            tau = qdot - M\(C*qdot + N)*self.timestep;
            z(1:numN) = Jn'*tau + (normals*q - alphas)/self.timestep;
            z(numN+1:numN+numT) = Jt'*tau;

            % Calculate the problem matrix for the QP
            P = zeros(numT+2*numN, numT+2*numN);
            u = zeros(numN,numT);
            
            for n = 1:numN
                u(n,numN*(n-1)+1:numN*n) = 1;
            end
            
            Mn = M\Jn;
            Mt = M\Jt;
           
            P(1:numN, 1:numN) = normals*Mn;
            P(1:numN, numN+1:numN+numT) = normals*Mt;
            P(numN+1:numN+numT, 1:numN) = Jt' * Mn;
            P(numN+1:numN+numT, numN+1:numN+numT) = Jt' * Mt;
            P(numN+1:numN+numT, numN+numT+1:end) = u';
            P(numN + numT+1:end,1:numN) = eye(numN) * self.friction;
            P(numN+numT+1:end,numN+1:numN+numT) = -u;
           
            % Solve the LCP
            [f,r] = self.solver(P,z);
            % Get the normal and tangential forces
            f = f(1:numN+numT);

        end
        function checkSolver(self)
            if strcmpi(self.contactSolver,'PATH')
                self.solver = @self.pathSolver;
            elseif strcmpi(self.contactSolver,'NCP')
                %self.opts.SpecifyObjectiveGradient = true;
                self.solver = @self.ncpSolver;
                self.func = @self.ncpFunc;
            elseif strcmpi(self.contactSolver,'SLCP')
                self.solver = @self.remSolver;
                %self.func = @self.slcpFunc;
            elseif strcmpi(self.contactSolver,'SLCP_Logistic')
                self.solver = @self.ncpSolver;
                self.func = @self.slcp_logistic;
            else
                error('Unknown Contact Solver');
            end
        end
        function [f, r] = ncpSolver(self,P,z)
            fun = @(x) self.func(x,P,z,self.sigma);
            x0 = 0*ones(size(z));
            %f = fminunc(fun,x0,self.opts);
            f = lsqnonlin(fun,x0,[],[],self.opts);
            %f = fmincon(fun,x0,-P,z,[],[],zeros(size(z)),[]);
            self.lastSoln = f;
            r = fun(f);
            r = sum(r.^2);
        end
        function [f,r] = remSolver(self,P,z)
           fun = @(x) self.slcpFunc(x,P,z,self.sigma);
           opt = optimoptions('fmincon','MaxFunctionEvaluations',5000,'StepTolerance',1e-8);
           if isempty(self.lastSoln)
               self.lastSoln = zeros(size(z));
           end
           %x0 = ones(size(z));
           %f = fminunc(fun,self.lastSoln,opt);
           f = fmincon(fun,self.lastSoln,-P,z,[],[],zeros(size(z)),[],[],opt);
           self.lastSoln = f;
           r = fun(f); 
        end
    end
      
    methods (Static)
        function [f, r] = pathSolver(P,z)
            try
                f = pathlcp(P,z);
                r = f'*(P*f + z);
            catch
                f = zeros(size(z));
                r = f'*(P*f + z);
                warning('PATH Solver failed');
            end
        end       
        function r = ncpFunc(x, P, z, ~)
           % Dummy Variables
           a = x;
           b = P*x+z;
           s = sqrt(a.^2 + b.^2);
           % Evaluate the NCP function
           phi = a + b - s;
           r = phi;
           % Calculate the residual
           %r = 1/2*sum(phi.^2);
           % Calculate the gradient
           %I = eye(numel(x));
           %dPhi = I + P - diag(1./s)*(diag(a) + diag(b)*P);
           %dr = dPhi'*phi;
        end
        function r = slcpFunc(x,A,z, sigma)
            a = x;
            b = A*x + z;
            % Calculate the PDF and CDF
            p = normpdf(a,b,sigma);
            P = normcdf(a,b,sigma);
            % Calculate the expected residual 
            r = a.^2 - sigma.^2 .* (a + b).*p + (sigma.^2 + b.^2 -a.^2).*P;
            % Sum-of-expected-squares
            %r = sqrt(r);
            r = sum(r);
        end
        function L = slcp_logistic(x, P, z, sigma)
           % Dummy variables
           a = x;
           b = P*x + z;
           p = exp((a-b)./sigma);
           % calculate the residuals
           L = a - sigma.*log(1 + p);
           % Return the sum of the squared residuals
%            r = 1/2*sum(L.^2);
%            % Calculate the gradient
%            I = eye(size(a));
%            dL = I -  diag(p./(1+p))*(I - P);
%            dr = dL'*L;
        end
    end
end