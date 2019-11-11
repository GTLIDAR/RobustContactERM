function [f,J,domerr] = mcp_adfunjac(z,jacflag)

global mcp_adfun;

if jacflag == 1
  [f,J] = evalJ(mcp_adfun, z);
else
  J = [];
  f = eval(mcp_adfun, z);
end;

domerr = 0;
return;

