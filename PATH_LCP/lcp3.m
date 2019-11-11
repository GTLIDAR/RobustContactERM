function [f,J,domerr]=lcp3(z,jacflag) 

% Input is z - the current point
%          jacflag = 1 is the jacobian is needed, 0 otherwise
%
% Return - f - the function evaluation at z
%          j - the jacobian evaluation
%          domerr - number of undefined functions at z (division by zero)

M = eye(3); q = [-1; 1; 0];

if jacflag == 1
  J = sparse(M);
else
  J = [];
end

f = M*z + q;
domerr = 0;

return;
