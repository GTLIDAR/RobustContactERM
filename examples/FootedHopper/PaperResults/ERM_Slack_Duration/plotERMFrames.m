function  plotERMFrames( )
%PLOTERMFRAMES Summary of this function goes here
%   Detailed explanation goes here

data = load('3Second_Combined_HeightSweep/ERM_HeightUncertainty_Sigma_5E-01/FootedHopper_Flat_3Second_Combined_HeightSweep.mat');


ref = data.guess.traj.x;
xref = ref.eval(ref.getBreaks());
xerm = data.soln.xtraj.eval(data.soln.xtraj.getBreaks());

refCOM = zeros(2,size(xref,2));
ermCOM = zeros(2,size(xref,2));

for n = 1:size(xref,2)
   refCOM(:,n) = data.plant.centerOfMass(xref(1:5,n));
   ermCOM(:,n) = data.plant.centerOfMass(xerm(1:5,n));
end

% Plot the frames
frameidx = [1,27, 45, 66, 101];
qref = xref(1:5,:);
qerm = xerm(1:5,:);
animationUtilities.visualizeFrames(data.plant, {qref, qerm}, frameidx, {'Reference','ERM'});

% overlay the COM trajectories
plot(refCOM(1,:),refCOM(2,:),'-','LineWidth',1.5,'DisplayName','RefCOM');
plot(ermCOM(1,:),ermCOM(2,:),'-', 'LineWidth',1.5,'DisplayName','ERMCOM');

end

