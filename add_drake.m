function add_drake()
%% ADD_DRAKE: adds drake files to the MATLAB search path
%
%   add_drake is a wrapper for the M-file in Drake 
%   "addpath_drake.m"
%
%   add_drake changes the current directory to the directory containing 
%   addpath_drake.m, calls the script, and then changes the directory
%   back to the directory from which this script was called

% Record the current working directory
present = pwd;

% Switch to the Drake directory to add Drake
drake_path = '/home/ldrnach/Projects/Drake/drake/drake';
cd(drake_path);
addpath_drake;

% Switch back to the working directorys
cd(present);

end