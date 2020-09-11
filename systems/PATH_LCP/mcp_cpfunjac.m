function [f,J,domerr] = mcp_cpfunjac(z,jacflag)

if jacflag
   [f,J,domerr] = mcpdf(z);
else
   J = [];
   [f,domerr] = mcpf(z);
end

return

