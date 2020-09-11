function [x, lambda] = lp_lcp(c,A,b,l,u,x,t)
% function [x, lambda] = lp_lcp(c,A,b,l,u,x,t)
%  

Big = 1e20;

if (nargin < 2)
  error('one argument required to lp_lcp(c)');
end

if (nargin < 2)
  A = [];
end

if (nargin < 3)
  b = [];
end

if (nargin < 4 | isempty(l))
  l = -Big*ones(length(c),1);
end
 
if (nargin < 5 | isempty(u))
  u =  Big*ones(length(c),1);
end

if (nargin < 6 | isempty(x))
  x = zeros(length(c),1);
end

if (nargin < 7 | isempty(t))
  t = -ones(length(b),1);
end

[x, lambda] = pathlcp(sparse(length(c),length(c)),c,l,u,x,A,b,t);

return;

