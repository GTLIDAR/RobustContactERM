[q,M] = lcp3(zeros(3,1),1);
x = pathlcp(M,q)
% soln is [1; 0; 0]
if ~ isempty(find(x - [1; 0; 0]))
  error ('bad solution x from pathlcp run 1');
end

x = pathlcp(speye(500),-ones(500,1));
if ~ isempty(find(x - ones(500,1)))
  error ('bad solution x from pathlcp run 2');
end

x = pathmcp(rand(3,1),zeros(3,1),inf*ones(3,1),'lcp3')
% soln is [1; 0; 0];
if ~ isempty(find(x - [1; 0; 0]))
  error ('bad solution x from pathmcp run 1');
end

x = pathmcp(rand(11,1),zeros(11,1),inf*ones(11,1),'mcp_funjac')
f = mcp_funjac(x,0)
resid = f'*x
if abs(resid) > 5e-6
  error ('bad solution x from pathmcp run 2');
end
display 'runtests completed OK'
