classdef UncertainDistanceERM < GaussianERM
    
   properties
        timestep;
        numN;
        numT;
   end
   methods
       function obj = UncertainDistanceERM()
           
       end
       function [m_mu, dmu_f, dmu_y] = ermMean(obj, f, P, w, dP, dw)
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
       end
       function [m_sigma, dsigma_f, dsigma_y]  = ermDeviation(obj, f, ~, ~, dP, ~)
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
           
           
           % Get the number of additional variables
           numY = size(dP, 3);

           % Calculate the variance
           m_sigma = obj.sigma ./ obj.timestep * ones(obj.numN,1);
           
           % The derivative of the variance with respect to the ERM
           % solution
           dsigma_f = zeros(obj.numN,obj.numN,length(f));
           % The derivative of the variance with respect to any other
           % parameters
           dsigma_y  = zeros(obj.numN,numY);
       end
   end
end