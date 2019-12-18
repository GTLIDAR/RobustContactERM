classdef UncertainFrictionERM < GaussianERM
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        numN;
        numT;
    end
    methods
        function obj = UncertainFrictionERM()
            
            
        end
        function [m_mu, dmu_f, dmu_y] = ermMean(obj, f, ~, ~, dP, ~)
            %% ERMMEAN: The mean of the Gaussian Distribution for the ERM problem
            %
            %   ermMean returns the mean of the Gaussian Distribution used
            %   in the ERM Contact problem
            %
            %   Arguments:
            %       OBJ:    an UncertainFrictionERM object
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
            
            % Get the number of additional variables
            numY = size(dP, 3);
            
            % Calculate the mean
            m_mu = obj.mu .* f(1:obj.numN,:) - e * f(obj.numN+1:obj.numN+obj.numT,:);
            % The derivative of the mean with respect to the ERM solution f
            dmu_f = [eye(obj.numN)*obj.mu, -e, zeros(obj.numN)];
            % The derivative with respect to any other parameters
            dmu_y = zeros(obj.numN, numY);
        end
        function [m_sigma, dsigma_f, dsigma_y]  = ermDeviation(obj, f, ~, ~, dP, ~)
            %% ERMDeviation: The standard deviation of the Gaussian Distribution for the ERM problem
            %
            %   ermDeviation returns the standard deviation of the Gaussian 
            %   Distribution underlying the ERM problem
            %
            %   Arguments:
            %       OBJ:    an UncertainFrictionERM object
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
            m_sigma = obj.sigma .* f(1:obj.numN,:);
            
            % The derivative of the variance with respect to the ERM
            % solution
            dsigma_f = zeros(obj.numN, length(f));
            dsigma_f(1:obj.numN, 1:obj.numN) = obj.sigma * eye(obj.numN); 
            % The derivative of the variance with respect to any other
            % parameters
            dsigma_y  = zeros(obj.numN, numY);
        end
    end
end

