function plotAllData(trajData,name, field, nominal)
% Helper function for plotting several trajectories at once
%
%   Luke Drnach
%   January 22, 2020


f1 = figure('name',[name,'- Trajectories']);  %Trajectories
f2 = figure('name',[name,'- Forces']);  %Contact forces
S = [1,0,0,0;0,1,-1,0];
% loop over all the data
for n = 1:length(trajData)
    % Get the points from the trajectory
    [t,x] = getPointsFromTrajectory(trajData(n).xtraj);
    [~,u] = getPointsFromTrajectory(trajData(n).utraj);
    [~,l] = getPointsFromTrajectory(trajData(n).ltraj);
    forces = S*l(:,1:end-1);
    u = u(:,1:end-1);
    if strcmpi(field, 'friction')
        lgdstr = sprintf('\\mu = %0.2f',trajData(n).friction);
    elseif strcmpi(field, 'sigma')
        lgdstr = sprintf('\\sigma = %0.4f',trajData(n).sigma);
    elseif strcmpi(field,'multiplier')
        lgdstr = sprintf('M = %0.2f',trajData(n).multiplier);
    else
        lgdstr = sprintf('Trial %d',n);
    end
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
    plot(t(1:end-1),u,'LineWidth',1.5);
    hold on;
    ylabel('Control');
    xlabel('Time');
    % Plot the reaction forces
    figure(f2);
    subplot(2,1,1);
    plot(t(1:end-1),forces(1,:),'LineWidth',1.5,'DisplayName',lgdstr);
    hold on;
    ylabel('Normal');
    title('Contact Forces');
    subplot(2,1,2);
    plot(t(1:end-1),forces(2,:),'LineWidth',1.5);
    hold on;
    ylabel('Tangential');
    xlabel('Time');
end

if nargin == 4
   % Plot the nominal trajectory
   [t, x] = getPointsFromTrajectory(nominal.xtraj);
   [~, u] = getPointsFromTrajectory(nominal.utraj);
   [~, f] = getPointsFromTrajectory(nominal.ltraj);
   forces = S*f(:,1:end-1);
   u = u(:,1:end-1);
    lgdstr = 'nominal';
        % Plot the state and control
    figure(f1);
    subplot(3,1,1);
    plot(t,x(1,:),'k-','LineWidth',1,'DisplayName',lgdstr);
    subplot(3,1,2);
    plot(t,x(3,:),'k-','LineWidth',1);
    subplot(3,1,3);
    plot(t(1:end-1),u,'k-','LineWidth',1);
    % Plot the reaction forces
    figure(f2);
    subplot(2,1,1);
    plot(t(1:end-1),forces(1,:),'k-','LineWidth',1,'DisplayName',lgdstr);
    subplot(2,1,2);
    plot(t(1:end-1),forces(2,:),'k-','LineWidth',1);
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