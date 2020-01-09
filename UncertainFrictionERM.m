classdef UncertainFrictionERM < GaussianERM
    %% UNCERTAINFRICTIONERM: Solves the ERM problem when there is Gaussian Uncertainty on the friction coefficient
    %
    %   UncertainFrictionERM is a concrete instance of GaussianERM that
    %   implements the mean and standard deviation of the slack variables
    %   of the LCP problem when the friction coefficient is uncertain and
    %   the uncertainty is normally distributed.
    
    %   Luke Drnach
    %   December 18, 2019
    
    properties
        numN;   %Number of normal force components
        numT;   %Number of tangential force components
    end
    methods
        function obj = UncertainFrictionERM(plant,mu,sigma)
            %% UncertainFrictionERM: Creates an instance of the UncertainFrictionERM class
            %
            %    UncertainFrictionERM is used when the shape of the terrain
            %    and the distance to the terrain is known, but the friction
            %    coefficient is uncertain, with the uncertainty modeled as
            %    a Gaussian distribution.
            %
            %   Arguments:
            %       PLANT:  An instance of a subclass of the
            %               DifferentiableContactDynamics class
            %       m:      Ns x 1 double, the constant means of the
            %               Gaussian distribution describing the parameter
            %               uncertainty.
            %       s:      Ns x 1 double, the standard deviations of the
            %               Gaussian distribution describing the parameter
            %               uncertainty.
            %               variables subject to uncertainty.
            %   Return Values:
            %       OBJ:     An instance of the GaussianERM class
            %
            %   Notes: Ns is the number of variables subject to uncertainty.
            %   In general, the number of TRUE values in IDX and Ns should
            %   be equal.
            
            % Initialize the parent class
            obj = obj@GaussianERM(mu,sigma,0);
            % Get the number of normal and tangential force components
            q = zeros(plant.numQ,1);
            [Jn,Jt] = plant.contactJacobian(q);
            obj.numN = size(Jn,1);
            obj.numT = size(Jt,1);
            % Record that the friction cone is uncertain
            obj.uncertainIdx = false(2*obj.numN + obj.numT,1);
            obj.uncertainIdx(obj.numN+obj.numT+1:end,:) = true;
        end
        function [m_mu, dmu_x, dmu_y, dmu_xx, dmu_xy] = ermMean(obj, x, P, ~, varargin)
            %% ERMMEAN: The mean of the Gaussian Distribution for the ERM problem
            %
            %   ermMean returns the mean of the Gaussian Distribution used
            %   in the ERM Contact problem
            %
            %   Arguments:
            %       OBJ:    an UncertainFrictionERM object
            %       x:      Nx1 double, the solution vector to the ERM
            %               problem
            %       P:      NxN double, the matrix in the LCP problem 
            %               x'(Px + w) = 0 
            %       w:      Nx1 double, the vector in the LCP problem
            %               x'(Px + w) = 0
            %       dP:     NxNxM double, an array of derivatives of P with
            %               respect to other variables (usually the state 
            %               and controls of a dynamical system) 
            %       dw:     NxM double, an array of derivatives of w with
            %       respect to other variables (usually the state and
            %               controls of a dynamical system)
            %   Return Values:
            %       m_mu:   numNx1 double, the mean of the Gaussian
            %               Distribution for the ERM problem. m_mu is the
            %               mean of the component of Px+w that encodes the
            %               friction cone constraint. numN is the number of
            %               normal friction components (the number of
            %               contact points).
            %       dmu_x:  numNxN double, the partial derivative of m_mu 
            %               with respect to the LCP solution, x
            %       dmu_y:  numNxM double, the partial derivative of m_mu 
            %               with respect to any other variables (usually
            %               state and controls of a dynamical system)
            %       dmu_xx: numNxNxN double, the second partial derivative
            %               of m_mu with respect to the LCP  solution x
            %       dmu_xy: numNxNxM double, the mixed partial derivatives
            %               of m_mu with respect to the LCP solution x and
            %               and the parameters y

            u = -P(obj.numN + obj.numT + 1:end, obj.numN + 1 : obj.numN + obj.numT);
            % Calculate the mean
            m_mu = obj.mu .* x(1:obj.numN,:) - u * x(obj.numN+1:obj.numN+obj.numT,:);
            % The derivative of the mean with respect to the ERM solution f
            dmu_x = [eye(obj.numN)*obj.mu, -u, zeros(obj.numN)];
            if nargout > 2
                dP = varargin{1};
                numY = size(dP,3);
                % The derivative with respect to any other parameters
                dmu_y = zeros(obj.numN, numY);
                % The first derivatives are constant, so the higher order
                % derivatives are zero
                dmu_xx = zeros(obj.numN, length(x), length(x));
                dmu_xy = zeros(obj.numN, length(x), numY);
            end
        end
        function [m_sigma, dsigma_x, dsigma_y, dsigma_xx, dsigma_xy]  = ermDeviation(obj, x, ~, ~, varargin)
            %% ERMDeviation: The standard deviation of the Gaussian Distribution for the ERM problem
            %
            %   ermDeviation returns the standard deviation of the Gaussian 
            %   Distribution underlying the ERM problem
            %
            %   Arguments:
            %       OBJ:    an UncertainFrictionERM object
            %       x:      Nx1 double, the solution vector to the ERM
            %               problem
            %       P:      NxN double, the matrix in the LCP problem 
            %               x'(Px + w) = 0 
            %       w:      Nx1 double, the vector in the LCP problem
            %               x'(Px + w) = 0
            %       dP:     NxNxM double, an array of derivatives of P with
            %               respect to other variables (usually the state 
            %               and controls of a dynamical system) 
            %       dw:     NxM double, an array of derivatives of w with
            %       respect to other variables (usually the state and
            %               controls of a dynamical system)
            %   Return Values:
            %       m_sigma:numNx1 double, the standard deviation of the Gaussian
            %               Distribution for the ERM problem. m_sigma is the
            %               standard deviation of the component of Pf+w that encodes the
            %               friction cone constraint. numN is the number of
            %               normal friction components (the number of
            %               contact points).
            %       dsigma_x:  numNxN double, the partial derivative of m_sigma 
            %               with respect to the LCP solution, x
            %       dsigma_y:  numNxM double, the partial derivative of m_sigma 
            %               with respect to any other variables (usually
            %               state and controls of a dynamical system)
            %       dsigma_xx: numNxNxN double, the second partial
            %               derivative of m_sigma with respect to the LCP
            %               solution x
            %       dsigma_xy: numNxNxM double, the mixed partial derivative
            %               of m_sigma with respect to the LCP solution x
            %               and the parameters y.
            numX = numel(x);
            % Calculate the variance
            m_sigma = obj.sigma .* x(1:obj.numN,:);
            % The derivative of the variance with respect to the ERM
            % solution
            dsigma_x = zeros(obj.numN, numX);
            dsigma_x(1:obj.numN, 1:obj.numN) = obj.sigma * eye(obj.numN); 
            if nargout > 2
                dP = varargin{1};
                numY = size(dP,3);
                % The derivative of the variance with respect to any other
                % parameters
                dsigma_y  = zeros(obj.numN, numY);
                % The first derivatives are constants, so all the higher order
                % derivatives are zero
                dsigma_xx = zeros(obj.numN, numX, numX);
                dsigma_xy = zeros(obj.numN, numX, numY);
            end
        end
    end
end

