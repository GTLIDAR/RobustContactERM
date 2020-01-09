classdef UncertainDistanceERM < GaussianERM
    %% UNCERTAINDISTNACEERM: Solves the ERM problem when there is Gaussian Uncertainty on the distance to the terrain
    %
    %  UncertainDistanceERM implements a concrete class of GaussianERM
    %  where the shape of the terrain is known but the distance to the
    %  terrain is uncertain and the uncertainty is normally distributed.
    
    %   Luke Drnach
    %   December 18, 2019
    
   properties
        timestep;   % Step size for the time-discretization of the contact dynamics
        numN;       % Number of normal force components
        numT;       % Number of tangential force components
   end
   methods
       function obj = UncertainDistanceERM(plant,mu,sigma)
           %% UncertainDistanceERM: Creates an instance of the UncertainDistanceERM class
           %
           %    UncertainDistanceERM is used when the shape of the terrain
           %    and the friction coefficient are known, but the distance
           %    to the terrain is uncertain, with the uncertainty be
           %    modeled as a Gaussian distribution.
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
           % Calculate the number of normal and tangential components
           q = zeros(plant.numQ,1);
           [Jn,Jt] = plant.contactJacobian(q);
           obj.numN = size(Jn,1);
           obj.numT = size(Jt,1);
           % Mark that the normal force component is uncertain
           obj.uncertainIdx = false(2*obj.numN + obj.numT,1);
           obj.uncertainIdx(1:obj.numN,:) = true;
           % Record the timestep for use later
           obj.timestep = plant.timestep;
       end
       function [m_mu, dmu_x, dmu_y, dmu_xx, dmu_xy] = ermMean(obj, x, P, w, varargin)
           %% ERMMEAN: The mean of the Gaussian Distribution for the ERM problem
           %
           %   ermMean returns the mean of the Gaussian Distribution used
           %   in the ERM Contact problem
           %
           %   Arguments:
           %       OBJ:    an UncertainDistanceERM object
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
           %               mean of the component of Pf+w that encodes the
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
           
           % Get the size of the problem and the number of additional
           % variables
           
           % Calculate the mean
           z = P*x + w;
           m_mu = z(1:obj.numN,:) - obj.mu ./ obj.timestep;
           
           % The derivative of the mean with respect to the ERM solution x
           dmu_x = zeros(obj.numN, numel(x));
           dmu_x(1:obj.numN, 1:obj.numN+obj.numT) = P(1:obj.numN, 1:obj.numN+obj.numT);
           
           if nargout > 2
               dP = varargin{1};
               dw = varargin{2};
               % The derivative with respect to any other parameters
               dP_x = times(dP, reshape(x, [1,numel(x), 1]));
               dP_x = sum(dP_x, 2);
               dP_x = reshape(dP_x, [size(dP_x,1),size(dP_x, 3)]);
               dz_y = dP_x + dw;
               dmu_y = dz_y(1:obj.numN,:);
               
               % Second derivative with respect to the ERM solution, f
               dmu_xx = zeros(obj.numN, numel(x), numel(x));
               
               % Mixed partial derivative with respect to the solution f and
               % the parameters y
               dmu_xy = dP(1:obj.numN, :,:);
           end
       end
       function [m_sigma, dsigma_x, dsigma_y, dsigma_xx, dsigma_xy]  = ermDeviation(obj, x, ~, ~, varargin)
           %% ERMDeviation: The standard deviation of the Gaussian Distribution for the ERM problem
           %
           %   ermDeviation returns the standard deviation of the Gaussian
           %   Distribution underlying the ERM problem
           %
           %   Arguments:
           %       OBJ:    an UncertainDistanceERM object
           %       x:      Nx1 double, the solution vector to the ERM
           %               problem
           %       P:      NxN double, the matrix in the LCP problem
           %               x('Px + w) = 0
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

           % Calculate the variance
           m_sigma = obj.sigma ./ obj.timestep * ones(obj.numN,1);
           
           % NOTE: Because m_sigma is constant, all the derivatives are
           % zero
           
           % The derivative of the variance with respect to the ERM
           % solution
           dsigma_x = zeros(obj.numN,length(x));
           if nargout > 2
               dP = varargin{1};
               numY = size(dP,3);
               numX = numel(x);
               % The derivative of the variance with respect to any other
               % parameters
               dsigma_y  = zeros(obj.numN, numY);
               % Second derivative
               dsigma_xx = zeros(obj.numN, numX, numX);
               % Mixed derivative
               dsigma_xy = zeros(obj.numN, numX, numY);
           end
       end
   end
end