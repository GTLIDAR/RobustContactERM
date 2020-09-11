function [x, lambda] = qp_lcp(Q,c,A,b,l,u,x,t)
% function [x, lambda] = qp_lcp(Q,c,A,b,l,u,x,t)
%  

Big = 1e20;

if (nargin < 2)
  error('two arguments required to qp_lcp(Q,c)');
end

if (nargin < 3)
  A = [];
end

if (nargin < 4)
  b = [];
end

if (nargin < 5 | isempty(l))
  l = -Big*ones(length(c),1);
end
 
if (nargin < 6 | isempty(u))
  u =  Big*ones(length(c),1);
end

if (nargin < 7 | isempty(x))
  x = zeros(length(c),1);
end

if (nargin < 8 | isempty(t))
  t = -ones(length(b),1);
end

[x, lambda] = pathlcp(Q,c,l,u,x,A,b,t);

return;

