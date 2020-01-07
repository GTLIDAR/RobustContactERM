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
       function [m_mu, dmu_f, dmu_y, dmu_ff, dmu_yf] = ermMean(obj, f, P, w, dP, dw)
           %% ERMMEAN: The mean of the Gaussian Distribution for the ERM problem
           %
           %   ermMean returns the mean of the Gaussian Distribution used
           %   in the ERM Contact problem
           %
           %   Arguments:
           %       OBJ:    an UncertainDistanceERM object
           %       f:      Nx1 double, the solution vector to the ERM
           %               problem
           %       P:      NxN double, the matrix in the LCP problem
           %               f'Pf + f'w = 0
           %       w:      Nx1 double, the vector in the LCP problem
           %               f'Pf+f'w = 0
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
           %       dmu_f:  numNxN double, the partial derivative of m_mu
           %               with respect to the LCP solution, f
           %       dmu_y:  numNxM double, the partial derivative of m_mu
           %               with respect to any other variables (usually
           %               state and controls of a dynamical system)
           %       dmu_ff: numNxNxN double, the second partial derivative 
           %               of m_mu with respect to the LCP  solution f
           %       dmu_fy: numNxNxM double, the mixed partial derivatives
           %               of m_mu with respect to the LCP solution f and 
           %               and the parameters y
           
           
           % Get the size of the problem and the number of additional
           % variables
           
           % Calculate the mean
           z = P*f + w;
           m_mu = z(1:obj.numN,:) - obj.mu ./ obj.timestep;
           
           % The derivative of the mean with respect to the ERM solution f
           dmu_f = zeros(obj.numN, numel(f));
           dmu_f(1:obj.numN, 1:obj.numN+obj.numT) = P(1:obj.numN, 1:obj.numN+obj.numT);
           % The derivative with respect to any other parameters
           dP_f = times(dP, reshape(f, [1,numel(f), 1]));
           dP_f = sum(dP_f, 2);
           dP_f = reshape(dP_f, [size(dP_f,1),size(dP_f, 3)]);
           dz_y = dP_f + dw;        
           dmu_y = dz_y(1:obj.numN,:);
           
           % Second derivative with respect to the ERM solution, f
           dmu_ff = zeros(obj.numN, numel(f), numel(f));
           
           % Mixed partial derivative with respect to the solution f and
           % the parameters y
           dmu_yf = dP(1:obj.numN, :,:);
       end
       function [m_sigma, dsigma_f, dsigma_y, dsigma_ff, dsigma_fy]  = ermDeviation(obj, f, ~, ~, dP, ~)
           %% ERMDeviation: The standard deviation of the Gaussian Distribution for the ERM problem
           %
           %   ermDeviation returns the standard deviation of the Gaussian
           %   Distribution underlying the ERM problem
           %
           %   Arguments:
           %       OBJ:    an UncertainDistanceERM object
           %       f:      Nx1 double, the solution vector to the ERM
           %               problem
           %       P:      NxN double, the matrix in the LCP problem
           %               f'Pf + f'w = 0
           %       w:      Nx1 double, the vector in the LCP problem
           %               f'Pf+f'w = 0
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
           %       dsigma_f:  numNxN double, the partial derivative of m_sigma
           %               with respect to the LCP solution, f
           %       dsigma_y:  numNxM double, the partial derivative of m_sigma
           %               with respect to any other variables (usually
           %               state and controls of a dynamical system)
           %       dsigma_ff: numNxNxN double, the second partial 
           %               derivative of m_sigma with respect to the LCP
           %               solution f
           %       dsigma_fy: numNxNxM double, the mixed partial derivative
           %               of m_sigma with respect to the LCP solution f
           %               and the parameters y.
           
           % Get the number of additional variables
           numY = size(dP, 3);

           % Calculate the variance
           m_sigma = obj.sigma ./ obj.timestep * ones(obj.numN,1);
           
           % NOTE: Because m_sigma is constant, all the derivatives are
           % zero
           
           % The derivative of the variance with respect to the ERM
           % solution
           dsigma_f = zeros(obj.numN,length(f));
           % The derivative of the variance with respect to any other
           % parameters
           dsigma_y  = zeros(obj.numN,numY);
           % Second derivative
           dsigma_ff = zeros(obj.numN, length(f), length(f));
           % Mixed derivative
           dsigma_fy = zeros(obj.numN, length(f), numY);
       end
   end
end