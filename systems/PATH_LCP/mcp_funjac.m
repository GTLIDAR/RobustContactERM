function [F, J, domerr] = mcp_funjac(z, jacflag)

% initialize
z = z(:);
F = [];
J = [];
domerr = 0;

% obtain variable values
x = z(1:6);     % quantities
p_d = z(7:9);   % demand price
p_s = z(10:11); % supply price

% check for domain violation
if (p_d <= 0)
  domerr = 1;
  return;
end

cost   = [ 0.225; 0.153; 0.162; 0.225; 0.162; 0.126; ];
demand = [ -325; -300; -275; ] ./ (p_d.^0.5);
supply = [ 325; 575; ];

A = sparse( [
   1   0   0   1   0   0;
   0   1   0   0   1   0;
   0   0   1   0   0   1;
  -1  -1  -1   0   0   0;
   0   0   0  -1  -1  -1;
] );

F = [ (-A'*[p_d; p_s] + cost); A*x + [demand; supply] ];

if (jacflag)
  J = [ sparse(6,6) -A'; 
        A diag( [ -0.5 * [-325; -300; -275;] ./ (p_d.^1.5); 0; 0] ) ];
end

return

