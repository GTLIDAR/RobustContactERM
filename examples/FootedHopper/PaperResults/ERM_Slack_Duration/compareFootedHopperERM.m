function compareFootedHopperERM(plant, ref, erm, labels, savename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Get the configuration trajectory and foot trajectory of the reference
nQ = plant.getNumPositions();
t = ref.xtraj.getBreaks();
x = ref.xtraj.eval(t);
qref = x(1:nQ,:);
footRef = zeros(2,length(t));
for k = 1:length(t)
    p = plant.kinematics(qref(:,k));
    footRef(:,k) = mean(p, 2);
end

%Configuration and foot trajectories of ERM
if ~iscell(erm)
    erm = {erm};
end
nERM = length(erm);
qerm = cell(1,nERM);
footERM = cell(1,nERM);
for n = 1:nERM
    x = erm{n}.xtraj.eval(t);
    qerm{n} = x(1:nQ,:);
    % Calculate foot positions
    footERM{n} = zeros(2,length(t));
    for k = 1:length(t)
        p = plant.kinematics(qerm{n}(:,k));
        footERM{n}(:,k) = mean(p,2);
    end
end

% Plot the base and foot trajectories
fig = figure();
subplot(2,1,1);
for n = 1:nERM
    plot(t,qerm{n}(2,:),'LineWidth',1.5,'DisplayName',labels{n});
    hold on;
end
plot(t,qref(2,:),'LineWidth',1.5,'DisplayName','Reference');
legend show;
legend boxoff;
ylabel('Base Height (m)');
subplot(2,1,2);
for n = 1:nERM
   plot(t,footERM{n}(2,:),'LineWidth',1.5);
   hold on;
end
plot(t,footRef(2,:),'LineWidth',1.5);
ylabel('Foot Height (m)');
xlabel('Time (s)');
savefig(fig, [savename,'_BaseFoot.fig']);
% Calculate the average change in foot height

avgFootRef = mean(footRef(2,:));
avgBaseRef = mean(qref(2,:));
FootInc = zeros(1,nERM);
BaseInc = zeros(1,nERM);
for n = 1:nERM
   FootInc(n) = 100*(mean(footERM{n}(2,:)) - avgFootRef) ./(avgFootRef);
   BaseInc(n) = 100*(mean(qerm{n}(2,:)) - avgBaseRef)./avgBaseRef;
   fprintf('%s mean foot clearance increase: %f%%\n', labels{n}, FootInc(n)); 
   fprintf('%s mean base clearance increase: %f%%\n', labels{n}, BaseInc(n));
end
fig = figure();
subplot(2,1,1);
p = plot(BaseInc,'o');
ylabel('% Increase Base Height')
p.MarkerFaceColor = p.Color;
set(gca,'XTick',1:nERM,'XTickLabels',labels,'XTickLabelRotation',20);
subplot(2,1,2);

p = plot(FootInc,'o');
ylabel('% Increase Foot Height');
p.MarkerFaceColor = p.Color;
set(gca,'XTick',1:nERM,'XTickLabels',labels,'XTickLabelRotation',20);
% Animate the trajectories
pause(1.0);
solns = {ref, erm{:}};
labels = ['Reference',labels];
animateMultiSolution(plant, solns, labels, [savename, '.avi']);
end

