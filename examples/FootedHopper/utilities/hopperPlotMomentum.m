function  hopperPlotMomentum(plant, soln)
%HOPPERPLOTMOMENTUM 
%   Plot the linear momentum, the change in linear momentum, and the
%   contact forces in one plot

% Luke Drnach
% August 20, 2020

t = soln.xtraj.getBreaks();
% Get the configuration and force trajectory
x = soln.xtraj.eval(t);
l = soln.ltraj.eval(t);
% Configuration
q = x(1:5,:);
% Velocity
dq = x(6:10,:);
% Contact forces
fN = [l(1,:);l(5,:)];
fT = [l(2,:)-l(3,:); l(6,:) - l(7,:)];

N = size(x,2);
p = zeros(2,N);

for n = 1:N
    M = plant.massMatrix(q(:,n));
    p(:,n) = M(1:2,:)*dq(:,n);
end
dp = [diff(p,1,2), zeros(2,1)];
h = diff(t);
dp(:,1:end-1) = dp(:,1:end-1);

% Gravity offset
g = 9.81*(plant.baseMass + sum(plant.masses));

figure();
ax = subplot(2,1,1);
yyaxis left;
plot(t, dp(1,:),'LineWidth',1.5);
ylabel('Momentum Change');
title('Horizontal Direction');

yyaxis right;
plot(t,sum(fT,1),'LineWidth',1.5);
ylabel('Friction Impulse');

matchAxes(ax);

ax = subplot(2,1,2);
yyaxis left;
plot(t,dp(2,:),'LineWidth',1.5);
ylabel('Momentum Change');
title('Vertical Direction');
yyaxis right;
plot(t,sum(fN,1),'LineWidth',1.5);
ylabel('Normal Impulse');
xlabel('Time (s)');

matchAxes(ax);

end

function ax = matchAxes(ax)
yyaxis(ax, 'left');
left_lims = ylim;
yyaxis(ax,'right');
right_lims = ylim;

ymin = min(left_lims(1), right_lims(1));
ymax = max(left_lims(2),right_lims(2));

yyaxis(ax, 'left');
ylim([ymin, ymax]);
yyaxis(ax,'right');
ylim([ymin, ymax]);
end

function ax = alignZero(ax)
yyaxis(ax, 'left');
left_lims = ylim;
yyaxis(ax,'right');
right_lims = ylim;

l_range = range(left_lims);
r_range = range(right_lims);

l_percent = left_lims/l_range;
r_percent = right_lims/r_range;

new_percent = [min(l_percent(1), r_percent(1)), max(l_percent(2), r_percent(2))];

yyaxis(ax, 'left')
ylim(new_percent * l_range);
yyaxis(ax, 'right');
ylim(new_percent * r_range);



end
