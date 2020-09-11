classdef RelaxedNonlinearComplementarityConstraint < CompositeConstraint
   % RelaxedNonlinearComplementarityConstraint:
   %    A constraint of the form z >= 0, f(x, z) >= 0, zf(x,z) = 0
   %
   % RelaxedNonlinearComplementarityConstraint implements all the modes of 
   %    NonlinearComplementarityConstraint_original, as well as a new mode
   %    for relaxing the constraint
   %
   %    mode 1: (default, standard constraint)
   %        z, f(x, z) >= 0
   %        z*f(x,z) = 0
   %
   %    mode 2: (slack variable method)
   %        z, gam >= 0
   %        gam - f(x, z) = 0
   %        z*gam = 0
   %
   %    mode 3: (FB Function)
   %        z + f - sqrt(z.^2 + f.^2) = 0        
   %
   %    mode 4: (prox function)
   %        z - max(0, z - r*f) = 0
   %
   %    mode 5: Relaxed constraint (new)
   %        s, z, f(x, z) >= 0
   %        s - zf(x,z) >= 0
   %
   % In the relaxed mode, s must be a decision variable, and the
   % corresponding objective min_s M*s must be added
   
   properties
      ncc_fun = []; 
      xdim;
      zdim;
   end
   
   methods
       function obj = RelaxedNonlinearComplementarityConstraint(fun, xdim, zdim, mode, slack)
          
           % Check the inputs
           if nargin < 4 || isempty(mode)
               mode = 1;
           end
           if nargin < 5 || isempty(slack)
               slack = 0;
           end          
           % Initialize using an empty constructor
           obj = obj@CompositeConstraint;
           obj.n_slack = 0;
           % Store the function handle
           obj.ncc_fun = fun;
           % Store the dimensions
           obj.xdim = xdim;
           obj.zdim = zdim;
           % Create the constraints
           switch mode
               case 1   
                   % Default
                   constraints{1} = BoundingBoxConstraint([-inf(xdim,1);zeros(zdim,1)],inf(zdim+xdim,1));
                   constraints{2} = FunctionHandleConstraint(zeros(zdim,1),inf(zdim,1),xdim+zdim,fun);
                   constraints{3} = FunctionHandleConstraint(zeros(zdim,1),zeros(zdim,1)+slack,xdim+zdim,@obj.prodfun);
               case 2
                   % Extra slack variables
                   obj.n_slack = zdim;
                   constraints{1} = BoundingBoxConstraint([-inf(xdim,1);zeros(2*zdim,1)],inf(2*zdim+xdim,1));
                   constraints{2} = FunctionHandleConstraint(zeros(zdim,1),zeros(zdim,1),xdim+2*zdim,@obj.slackeq);
                   constraints{3} = FunctionHandleConstraint(zeros(zdim,1),zeros(zdim,1)+slack,xdim+2*zdim,@obj.slackprod);
                   constraints{1} = constraints{1}.setName(sprintf('NonlinearLCPBoundingBoxConstraint-s'));
                   constraints{2} = constraints{2}.setName(sprintf('NonlinearLCPSlackEq-s'));
                   constraints{3} = constraints{3}.setName(sprintf('NonlinearLCPSlackProd-s'));
               case 3
                   % FB Function
                   constraints = FunctionHandleConstraint(zeros(zdim,1),zeros(zdim,1),xdim+zdim,@obj.fbfun);
               case 4
                   % Prox/MIN function
                   constraints = FunctionHandleConstraint(zeros(zdim,1),zeros(zdim,1),xdim+zdim,@obj.proxfun);
               case 5
                   % Relaxed constraint
                   constraints{1} = BoundingBoxConstraint([-inf(xdim, 1); zeros(zdim+1, 1)], inf(xdim + zdim + 1, 1));
                   constraints{2} = FunctionHandleConstraint(zeros(zdim,1), inf(zdim,1), xdim + zdim + 1, @obj.relaxfun);
                   constraints{3} = FunctionHandleConstraint(zeros(zdim,1), inf(zdim,1), xdim + zdim + 1, @obj.relaxprod);
                   constraints{1} = constraints{1}.setName(sprintf('NCC_BoundingBoxConstraint'));
                   constraints{2} = constraints{2}.setName(sprintf('NCC_FunctionConstraint'));
                   constraints{3} = constraints{3}.setName(sprintf('NCC_RelaxedProductConstraint'));
           end
           % Add the constriants to the object
           obj = obj.addConstraints(constraints);

       end
   end
   methods (Access = protected)
       function [f, df] = prodfun(obj, y)
           % Get the variables and the NCC Function value
           z = y(obj.xdim+1:obj.xdim + obj.zdim);
           [g, dg] = obj.ncc_fun(y);
           % The product - ELEMENTWISE
           f = z.*g;
           % The derivative
           df = diag(z)*dg + [zeros(obj.zdim, obj.xdim), diag(g)];
       end
       function [f, df] = slackeq(obj, y)
           % Get the variables
           x = y(1:obj.xdim);
           z = y(obj.xdim+1:obj.xdim+obj.zdim);
           s = y(obj.xdim+obj.zdim+1:end);
           % Get the function value
           [f, df] = obj.ncc_fun([x;z]);
           % Now apply the slack variable
           f = f - s;
           % The derivative
           df = [df zeros(obj.zdim)] - [zeros(obj.zdim, obj.zdim+obj.xdim) eye(obj.zdim)];
       end
       function [f, df] = slackprod(obj, y)
           z = y(obj.xdim+1:obj.xdim+obj.zdim);
           s = y(obj.xdim+obj.zdim+1:end);
           % The product of the decision and slack variables
           f = z .* s;
           % The derivative
           df = [zeros(obj.zdim , obj.xdim), diag(s), diag(z)];
       end
       function [f,df] = fbfun(obj, y)
           %Fisher-Burmeister function implementation
           x = y(1:obj.xdim);
           z = y(obj.xdim+1:obj.xdim+obj.zdim);
           % Evaluate the underlying function
           [g,dg] = obj.ncc_fun([x;z]);
           % The NCP function 
           f = z + g - sqrt(z.^2 + g.^2);
           % The derivative
           df = [zeros(obj.zdim, obj.xdim), eye(obj.zdim)] + dg - diag(1./(sqrt(z.^2 + g.^2 + 1e-6))) * ([zeros(obj.zdim, obj.xdim), diag(z)] + diag(g) * dg);           
       end
       function [f, df] = proxfun(obj, y)
          % The MIN (Max) NCP function implementation
          x = y(1:obj.xdim);
          z = y(obj.xdim + 1: obj.xdim+obj.zdim);
          % Evaluate the underlying function
          [g, dg] = obj.ncc_fun([x;z]);
          % The NCP function
          r = 1;
          f = z - max(0, z - r * g);
          % The derivative - assuming z < g
          df = [zeros(obj.zdim, obj.xdim), eye(obj.zdim)];
          % indicator function
          I_pos = find(z - r*g >=0);
          % The derivative for z > g
          df(I_pos, obj.zdim + I_pos) = 0;
          df(I_pos, :) = df(I_pos, :) - r * dg(I_pos, :);
       end
       function [f, df] = relaxfun(obj, y)
           x = y(1:obj.xdim);
           z = y(obj.xdim+1:obj.xdim+obj.zdim);
           
           [g, dg] = obj.ncc_fun([x;z]);
           f = g;
           df = [dg, zeros(obj.zdim, 1)];
       end
       function [f, df] = relaxprod(obj,y)
           % Get the variables
           x = y(1:obj.xdim);          
           z = y(obj.xdim+1:obj.xdim+obj.zdim);
           s = y(obj.xdim+obj.zdim+1);    % Slack
           % NCC Function evaluation
           [g, dg] = obj.ncc_fun([x;z]);
           % The product
           f = s - z.*g;
           % Derivative (df/dx, df/dz, df/ds)
           df = [-z.*dg, ones(obj.zdim,1)] - [zeros(obj.zdim, obj.xdim), diag(g), zeros(obj.zdim, 1)];
       end
%        function [f, df] = relaxCost(obj, y)
%            % Get the slack variable only
%            f = obj.mult * y(obj.xdim + obj.zdim + 1);
%            % Derivative
%            df = [zeros(1, obj.xdim + obj.zdim), obj.mult];
%        end
   end
end