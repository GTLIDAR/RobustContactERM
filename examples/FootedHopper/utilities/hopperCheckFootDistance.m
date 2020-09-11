function hopperCheckFootDistance(plant, soln)
%HOPPERCHECKFOOTDISTANCE Summary of this function goes here
%   Detailed explanation goes here

%   July 16, 2020
%   Luke Drnach


t = soln.xtraj.getBreaks();
x = soln.xtraj.eval(t);
q = x(1:5,:);

pToe = zeros(2,numel(t));
pHeel = zeros(2,numel(t));

for n = 1:numel(t)
   p = plant.kinematics(q(:,n));
   pToe(:,n) = p(:,1);
   pHeel(:,n) = p(:,2);
end

% Draw the foot trajectory
figure();
subplot(2,1,1);
plot(pToe(1,:),pToe(2,:),'LineWidth',1.5,'DisplayName','Toe');
hold on;
plot(pHeel(1,:),pHeel(2,:),'LineWidth',1.5,'DisplayName','Heel');

% Draw the terrain
xl = xlim();
[xterrain,yterrain] = plant.terrain.draw(xl,100);
plot(xterrain, yterrain, 'k', 'LineWidth',1.5,'DisplayName','Terrain');
legend show
legend boxoff

% Check the Normal Distance
phi = zeros(2,numel(t));
for n = 1:numel(t)
   phi(:,n) = plant.contactConstraints(q(:,n)); 
end
subplot(2,1,2);
plot(t,phi(1,:),'LineWidth',1.5);
hold on;
plot(t,phi(2,:),'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Normal Distance');

end

