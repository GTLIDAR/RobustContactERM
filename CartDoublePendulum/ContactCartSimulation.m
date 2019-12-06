function ContactCartSimulation()
%DIFFERENTIABLECONTACTCARTSIMULATION Summary of this function goes here
%   Detailed explanation goes here

% Add the necessary paths
here = pwd;
cd ..
add_drake;
cd(here);
addpath(genpath('PATH_LCP'));

name = 'ContactCart_Sim';

plant = ContactDrivenCart();
plant.cartHeight = 1.5;
plant.timestep = 0.0025;
plant.terrain.friction_coeff = 0.5;

x0 = [0, pi/4, 0, 0, 0, 0]';

Tf = 5;

[~, x] = plant.simulate(Tf, x0);


figure();
plot(0,0);
ax = gca;

draw = @(ax, x) plant.draw(x, ax);
utilities.animator(ax, draw, x(1:3,:),[name,'.avi']);
end

