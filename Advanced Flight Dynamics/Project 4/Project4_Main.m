clc; clear; close all;
clear Riccati; % Recalculate cached gains for each simulation run.
render_orbits = true; % Open the interactive Earth/formation scene after simulation.
%%% initial conditions %%%

rho_i = [10000;20000;2000]; %initial relative position, m
drho_i = [0.9;-1;0.5]; %initial relative velocity, m/s

rho_f = [50;50;0]; %desired final relative positio, m
drho_f = [0;0;0]; %desired final relative velocity m/s

yd = [rho_f;drho_f];
mf = 10; %follower mass in kg


tf = ceil(2*pi/sqrt(3.986e14/(6898e3)^3)); % One leader orbit (~5702 s); use 1000 for a short run.
load_system('Project4_Sim.slx');
set_param('Project4_Sim/Integrator', 'InitialCondition', '[rho_i;drho_i]');
result = sim('Project4_Sim.slx', 'ReturnWorkspaceOutputs', 'on');
%%
time = result.x.Time;
x = result.x.Data;
y = result.y.Data;
U = mf*result.U.Data; %changed to mass

x_rho = x(:,1:3); x_drho = x(:,4:6);

y_rho = y(:,1:3); % Output contains the three position components.
%%

% States Plot
figure(1);
subplot(3,1,1); plot(time,x_rho,'LineWidth',0.8); hold on;
xline(500); hold on; yline(5+50); hold on; yline(-5+50);
legend('x','y','z'); title('a) Relative Position');
xlabel('Time (s)'); ylabel('Position (m)');

subplot(3,1,2); plot(time,x_drho,'LineWidth',0.8); hold on; xline(500);
legend('x','y','z'); title('b) Relative Velocity');
xlabel('Time (s)'); ylabel('Velocity (m/s)');

subplot(3,1,3); plot(time,U,'LineWidth',0.8); hold on; 
xline(500); 
legend('x','y','z'); title('c) Control');
xlabel('Time (s)'); ylabel('Force (N)');

%% Render the leader and follower using the simulated LVLH states.
if render_orbits
    Render_Orbits(time, x, rho_f);
end







