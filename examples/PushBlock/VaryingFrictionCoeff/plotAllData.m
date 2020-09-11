function plotAllData(trajData,name)
% Helper function for plotting several trajectories at once
%
%   Luke Drnach
%   January 22, 2020


f1 = figure('name',[name,'- Trajectories']);  %Trajectories
f2 = figure('name',[name,'- Forces']);  %Contact forces

% loop over all the data
for n = 1:length(trajData)
    % Get the points from the trajectory
    [t,x] = getPointsFromTrajectory(trajData(n).xtraj);
    [~,u] = getPointsFromTrajectory(trajData(n).utraj);
    lgdstr = sprintf('\\mu = %0.2f',trajData(n).friction);
    % Plot the state and control
    figure(f1);
    subplot(3,1,1);
    plot(t,x(1,:),'LineWidth',1.5,'DisplayName',lgdstr);
    hold on;
    ylabel('Position');
    subplot(3,1,2);
    plot(t,x(3,:),'LineWidth',1.5);
    hold on;
    ylabel('Velocity');
    subplot(3,1,3);
    plot(t,u,'LineWidth',1.5);
    hold on;
    ylabel('Control');
    xlabel('Time');
    % Plot the reaction forces
    figure(f2);
    subplot(2,1,1);
    plot(t,trajData(n).forces(1,:),'LineWidth',1.5,'DisplayName',lgdstr);
    hold on;
    ylabel('Normal');
    title('Contact Forces');
    subplot(2,1,2);
    plot(t,trajData(n).forces(2,:),'LineWidth',1.5);
    hold on;
    ylabel('Tangential');
    xlabel('Time');
end
figure(f1);
subplot(3,1,1);
legend show;
legend boxoff;
figure(f2);
subplot(2,1,1);
legend show;
legend boxoff;

end