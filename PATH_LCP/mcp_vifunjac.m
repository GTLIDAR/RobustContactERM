function [f,J,domerr]=mcp_vifunc(z,jacflag) 

global mcp_vifunc;
global mcp_viconn;
global mcp_viconm;
global mcp_viconA;
global mcp_viconb;

if jacflag == 1
  [f,J,domerr] = feval(mcp_vifunc,z(1:mcp_viconn),jacflag);
  if (mcp_viconm > 0)
    f = [f - mcp_viconA'*z(mcp_viconn+1:mcp_viconn+mcp_viconm); 
         mcp_viconA*z(1:mcp_viconn) - mcp_viconb];
    J = [J, -mcp_viconA'; mcp_viconA, sparse(mcp_viconm,mcp_viconm)];
  end
else
  [f,J,domerr] = feval(mcp_vifunc,z(1:mcp_viconn),jacflag);
  f = [f - mcp_viconA'*z(mcp_viconn+1:mcp_viconn+mcp_viconm); 
       mcp_viconA*z(1:mcp_viconn) - mcp_viconb];
end

return;
