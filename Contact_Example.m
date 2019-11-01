% Set up the falling rod example
model = FallingRod(1, 0.5, 0.05, 0.002, 0.6, 0.0025);
X0 = [0.0, 1.0, pi/6, 0.0, 0.0, 4.0]';

% Set the solver option
model.contactSolver = 'PATH';
% Forward simulate a trajectory
[t, x, f, r] = model.simulate(X0, 1.0);

% Display the results
model.plotTrajectory(t,x,f);
frames = model.animate(x);

% w = VideoWriter('fallingrod.avi');
% open(w);
% for n = 1:length(frames)
%     writeVideo(w,frames(n))
% end
% close(w);