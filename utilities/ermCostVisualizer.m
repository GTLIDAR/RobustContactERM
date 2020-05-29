function ermCostVisualizer(plant, q, t, x, mu, sigma, x_label, mu_label)
%ERMCOSTVISUALIZER Summary of this function goes here
%   Detailed explanation goes here

if nargin < 7
   x_label = 'x'; 
   mu_label = '\mu';
end
% Open a video writer
writer = VideoWriter('ermCostVisualization.avi');
writer.open();
% Calculate the ERM Cost trace:
f_trace = ermCostGaussian(x, mu, sigma);
% Determine the size of the ERM window
xMax = max(abs(x));
mMax = max(abs(mu));
xs = linspace(-xMax, xMax, 1000);
mus = linspace(-mMax, mMax, 1000);
% Make a meshgrid for making the ERM Cost Map
[Xg, Mg] = meshgrid(xs, mus);
F_map = ermCostGaussian(Xg, Mg, sigma(1)*ones(size(Xg)));
% Make the initial plots
figure('units','normalized','outerposition',[0,0,1,1]);
% The ERM Costmap plot
subplot(3,2,[1,3]);
img = imagesc(xs, mus, log10(F_map));
axis xy;
xlabel(x_label);
ylabel(mu_label);
title('log10(ERM)');
colormap jet;
colorbar;
hold on;
img_marker = plot(x(1), mu(1), 'ko','MarkerFaceColor','k');

% The X plot
subplot(3,2,2);
plot(t, x, 'LineWidth',1.5);
ylabel(x_label);
xlim([t(1),t(end)]);
hold on;
fn_marker = plot(t(1), x(1),'ko','MarkerFaceColor','k');

% The MU plot
subplot(3,2,4);
plot(t, mu, 'LineWidth',1.5);
ylabel(mu_label);
xlim([t(1),t(end)]);
hold on;
dis_marker = plot(t(1), mu(1), 'ko','MarkerFaceColor','k');

% The ERM Cost plot
subplot(3,2,6);
plot(t, f_trace, 'LineWidth',1.5);
xlim([t(1),t(end)]);
ylabel('ERM Cost');
hold on;
erm_marker = plot(t(1), f_trace(1), 'ko','MarkerFaceColor','k');

% The system configuration plot
ax = subplot(3,2,5);
[config_frames, xlims, ylims] = animationframes(plant,q);
% Adapt the x-axis to make it 'square'
ax.Units = 'pixels';
pos = ax.Position;
r = pos(4)/pos(3);
yr = ylims(2) - ylims(1);
ylims(2) = ylims(1) + yr*r;
% Draw the terrain and the first configuration
[xterrain, yterrain] = plant.terrain.draw(xlims, size(q,2));
plot(xterrain, yterrain, 'k-', 'LineWidth',2);
hold on;
% Draw the first point on the trajectory
sketch = plot(config_frames{1}(1,:), config_frames{1}(2,:),'b-','LineWidth',2);
% Set the axis limits
set(ax, 'XLimMode','manual','YLimMode','manual','xlim',xlims,'ylim',ylims);

% Update loop!
for n = 1:length(x)
   % Update the image map
    F_map = ermCostGaussian(Xg, Mg,sigma(n)*ones(size(Xg)));
    set(img, 'CData', log10(F_map));
    % Update the marker positions
    set(img_marker,'XData',x(n),'YData',mu(n));
    set(fn_marker,'XData',t(n),'YData',x(n));
    set(dis_marker,'XData',t(n),'YData',mu(n));
    set(erm_marker,'XData',t(n),'YData',f_trace(n));
    % Update the animation frame
    set(sketch,'XData',config_frames{n}(1,:),'YData',config_frames{n}(2,:));
    set(ax, 'XLimMode','manual','YLimMode','manual','xlim',xlims,'ylim',ylims);
    drawnow;
    pause(0.1);
    % Get the current frame
    frame = getframe(gcf);
    % Write the current frame
    writer.writeVideo(frame);
end
writer.close();
end

function f = ermCostGaussian(x, mu, sigma)

% Calculate the pdf values
pdf = normpdf(x, mu, sigma);
cdf = normcdf(x, mu, sigma);
% Set the values when sigma == 0
pdf(sigma == 0) = 0;
cdf(and(sigma == 0, x <= mu)) = 0;
cdf(and(sigma == 0, x > mu)) = 1;

% Calculate the ERM function evaluation
f = x.^2 - sigma.^2 .* (x + mu) .* pdf + (sigma.^2 + mu.^2 - x.^2).*cdf;

end


function [frames, xlims, ylims] = animationframes(model, seq)

% Draw all the components of the trajectory first
T = size(seq, 2);
frames = cell(1,T);
xlims = zeros(1,2);
ylims = zeros(1,2);
for n = 1:T
    [x,y] = model.draw(seq(:,n));
    frames{n} = [x;y];
    xlims(1) = min([xlims(1), x]);
    xlims(2) = max([xlims(2),x]);
    ylims(1) = min([ylims(1),y]);
    ylims(2) = max([ylims(2),y]);
end
% Extend the limits by 5%
xlims(1) = xlims(1) * (1 - sign(xlims(1)) * 0.05);
xlims(2) = xlims(2) * (1 + sign(xlims(2)) * 0.05);
ylims(1) = ylims(1) * (1 - sign(ylims(1)) * 0.05);
ylims(2) = ylims(2) * (1 + sign(ylims(2)) * 0.05);
% Make the limits square
xr = range(xlims);
yr = range(ylims);
if xr > yr
    ylims = ylims * (xr / yr);
else
    xlims = xlims * (yr/xr);
end

end