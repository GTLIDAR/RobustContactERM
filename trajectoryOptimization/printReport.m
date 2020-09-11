function printReport(plant, soln, solverOptions, optimOptions, guesstag, filename)

if nargin < 6
   fileID = fopen('report.txt','w');
else
   fileID = fopen(filename,'w'); 
end
% Write the date and author line
fprintf(fileID, 'Date: %s\n', date);
fprintf(fileID, 'Author: %s\n',char(java.lang.System.getProperty('user.name')));
fprintf(fileID, '\n');

% Print out the class file
fprintf(fileID, 'Plant model: %s\n', class(plant));

%Print out the optimization parameters
fprintf(fileID, '\nRun Settings:\n');
fprintf(fileID, '\tDuration: [%0.2f, %0.2f]\n', optimOptions.duration(1), optimOptions.duration(2));
fprintf(fileID, '\tNumber of Knot Points: %d\n', optimOptions.nPoints);
fprintf(fileID, '\tIncluded Boundary Conditions? %s\n', string(isfield(optimOptions, 'stateConstraints')));
fprintf(fileID, '\tIncluded Running Cost? %s\n', string(isfield(optimOptions, 'runningCost')));
fprintf(fileID, '\tIncluded Final Cost? %s\n', string(isfield(optimOptions, 'finalCost')));

if nargin == 5
   fprintf(fileID, '\tTrajectory Initialization: %s\n', guesstag); 
end

% Print out the solver options
fprintf(fileID, '\nSNOPT Settings:\n');

fields = fieldnames(solverOptions);
for n = 1:length(fields)
    fprintf(fileID, '\t %s: %.2e\n', fields{n}, solverOptions.(fields{n}));
end
% Print out the optimizer options
fprintf(fileID, '\nRobustContactImplicitTrajectoryOptimizer Settings: \n');
fields = fieldnames(optimOptions.options);
for n = 1:length(fields)
   fprintf(fileID, '\t %s: %d\n', fields{n}, optimOptions.options.(fields{n})); 
end

% Print out the results of the optimization
fprintf(fileID, '\nSNOPT terminated after %4.2f seconds with exit status %d\n', soln.elapsed, soln.info);

if ~isempty(soln.infeasible)
   infeasible = unique(soln.infeasible);
   fprintf(fileID,'\nInfeasible Constraints:\n');
   for n = 1:length(infeasible)
      fprintf(fileID,'\t%s\n',infeasible{n}); 
   end
end

if optimOptions.options.nlcc_mode == 5
   fprintf(fileID,'\nMaximum NCC Relaxation: %8.2e\n',max(soln.slacks)); 
end

fprintf(fileID, '\nNotes:\n');
% Close the file
fclose(fileID);
end