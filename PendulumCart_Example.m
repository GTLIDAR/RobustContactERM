%% PendulumCart Example

name = 'PendulumCartExample_SLCP_Logistic';

% Create the pendulumDrivenCart model
model = PendulumDrivenCart();
model.contactSolver = 'SLCP_Logistic';
model.sigma = 1;
model.initialHeight = 1.5;
model.timestep = 0.01;
model.friction = 0.5;
% Set an initial condition
x0 = [0;pi/4;0;0;0;0];
model.draw(x0(1:3));


check = input('Continue (Y/N)?','s');
if strcmpi(check,'y')
    % Run the simulation
    T = 10;      % One second simulation
    [t,x,f,r] = model.simulate(x0,T);
    ax= gca;
    % Animate the simulation    
    draw = @(ax,x) model.draw(x,ax);
    %utilities.animator(ax,draw,x(1:3,:));
    utilities.animator(ax,draw,x(1:3,:),[name,'.avi']);
    % Make a figure for the diagnostics
    utilities.diagnostics(model, t,x(1:3,:),f(1,:),name);
    % Make a figure for the residuals
    utilities.errorPlots(t,f,r,name);
end

