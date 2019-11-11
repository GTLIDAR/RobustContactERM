function [y] = semimcplib(modelname)

[y, ylb, yub] = mcpinit(modelname);

[y] = semimcp(y,ylb,yub,'mcp_cpfunjac');

mpecclose;

return;

