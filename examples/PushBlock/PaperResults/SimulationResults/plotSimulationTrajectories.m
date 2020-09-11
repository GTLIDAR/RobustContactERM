function plotSimulationTrajectories()
%PLOTSIMULATIONTRAJECTORIES Summary of this function goes here
%   Detailed explanation goes here

%   Luke Drnach
%   February 10, 2020

% Load the simulation results
results = load('openLoopSimResults.mat');
results = results.sim_results;
% Concatentate the source and uncertainty together into one 'condition'
% string
conditions = cell(1,length(results));
rowLabel = cell(1,length(results));
for n = 1:length(results)
   if strcmpi(results(n).source, 'ERM')
      conditions{n} = sprintf('ERM \\sigma = %0.4f',results(n).uncertainty); 
      rowLabel{n} = sprintf('ERM_%0.4f', results(n).uncertainty);
   else
      conditions{n} = results(n).source; 
      rowLabel{n} = results(n).source;
   end
end
% Get the unique conditions
unique_conditions = unique(conditions,'stable');
rowLabel = unique(rowLabel,'stable');
%Sigma = Sigma(id);
% Initialize error arrays
xidx = 1:length(unique_conditions);
MeanError = zeros(numel(xidx),1);
MaxError = MeanError;
MinError = MeanError;
Number = MeanError;
% Plot the simulated position and velocity trajectories
for n = 1:length(unique_conditions)
   figure(); 
   idx = find(strcmpi(conditions, unique_conditions{n})); 
   for k = 1:length(idx)
        subplot(2,1,1);
        plot(results(idx(k)).time, results(idx(k)).state(1,:),'Linewidth',1.5,'DisplayName',sprintf('\\mu = %0.2f',results(idx(k)).friction));
        hold on;
        subplot(2,1,2);
        plot(results(idx(k)).time, results(idx(k)).state(3,:),'Linewidth',1.5);
        hold on;
   end
   % Add in the target position and labels
   subplot(2,1,1);
   plot([0, results(idx(k)).time(end)], results(idx(k)).target(1)*[1,1],'k-','DisplayName','target');
   ylabel('Position');
   title(unique_conditions{n});
   legend('show');
   legend('location','best');
   legend('boxoff');
   subplot(2,1,2);
   plot([0, results(idx(k)).time(end)], results(idx(k)).target(3)*[1,1],'k-');
   ylabel('Velocity');
   xlabel('Time (s)');
   
   % Plot the mean and range of the terminal error
   err = [results(idx).terminalError];
   err = err(1,:);
   MeanError(n) = mean(err);
   MaxError(n) = max(err);
   MinError(n) = min(err);
   Number(n) = numel(idx);
end
% Plot the mean and range of the terminal error in each condition
figure();
errorbar(xidx, MeanError, MaxError - MeanError, MinError - MeanError, 'ko', 'MarkerFaceColor','k','LineWidth',1.5);
ylabel('Final Position Error');
set(gca,'XTick',xidx, 'XTickLabel',unique_conditions, 'XTickLabelRotation',30);
xlim([0,n+1]);
% Make and store a table of all the position errors

errorTable = table(MeanError, MaxError, MinError,  Number, 'RowNames',rowLabel);
disp(errorTable);
save('TerminalPositionErrors.mat','errorTable');
end

