function terrainCheck(terrain)
%TERRAINCHECK Summary of this function goes here
%   Detailed explanation goes here




% Create a high-fidelity drawing of the terrain
[x_t, y_t] = terrain.draw(terrain.step_location + [-terrain.step_width,0], 100);
[x_s, y_s] = terrain.draw(terrain.step_location + [0, terrain.step_width], 100);
[x_h, y_h] = terrain.draw(terrain.step_location + terrain.step_width*[1,2], 100);

x = [x_t(:); x_s(:); x_h(:)];
y = [y_t(:); y_s(:); y_h(:)];

figure();
plot(x,y,'k', 'LineWidth',1.5);

yl = ylim;
yr = range(yl);
yl = yl + 0.1*[-yr, yr];

ylim(yl);

% Now create a set of sample points
yl = ylim;
xl = xlim;

% [xgrid, ygrid] = meshgrid(linspace(xl(1), xl(2), 10), linspace(yl(1), yl(2), 10));
 hold on;
% 
% dots = plot(xgrid(:), ygrid(:), 'bo');
% dots.MarkerFaceColor = dots.Color;
% dots.MarkerSize = 3;

% Add 'fake points' for the nearest point and the distance circle
x = mean(xl);
y = mean(yl);
plot(x,y,'o');
plot(x,y,'d');
plot(x,y,'-');

% Add the callback for clicking on the figure
set(gca,'XLimMode','manual','YLimMode','manual');
set(gca,'ButtonDownFcn',{@clickCallback, terrain});
set(gca,'ButtonDownFcn',{@clickCallback, terrain});

end

function clickCallback(obj, event, terrain)


point = event.IntersectionPoint(1:2);
plt = get(obj,'Children');
drawNearest(terrain, point, plt);
end

function plt = drawNearest(terrain, point, plt)
% Get and draw the nearest point solution

xN = terrain.nearest(point);
set(plt(3),'XData',point(1),'YData',point(2));
set(plt(2),'XData',xN(1),'YData',xN(2));
% Draw the circle around the current point
r = norm(point' - xN);
xc = r*cos(0:0.01:2*pi) + point(1);
yc = r*sin(0:0.01:2*pi) + point(2);
set(plt(1), 'XData', xc, 'YData', yc);
drawnow;
end


