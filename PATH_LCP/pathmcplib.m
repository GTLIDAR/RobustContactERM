function [y] = pathmcplib(modelname)

[y, ylb, yub] = mcpinit(modelname);

[y] = pathmcp(y,ylb,yub,'mcp_cpfunjac');

mcpclose;

return;

