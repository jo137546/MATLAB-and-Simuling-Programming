% Project 3 - Option A only: nonlinear two-satellite relative motion.
% Run this script with all supplied files in the same folder.
% Units throughout: km, seconds, radians (convert degrees only at input).
clear; clc; close all;
project_dir = fileparts(mfilename('fullpath'));
addpath(project_dir);

%% Constant parameters and assignment initial conditions
mu = 3.986e5; % Earth gravitational parameter, km^3/s^2 (Oct. 1 notes)
a_L = 8000; e_L = 0.05;
i_L = deg2rad(45); bigW_L = deg2rad(20); wumg_L = deg2rad(90);
nu_L0 = deg2rad(-10);
a_F = 8000; e_F = 0.05;
i_F = deg2rad(45); bigW_F = deg2rad(20); wumg_F = deg2rad(90);
nu_F0 = deg2rad(10);

orb_period = 2*pi*sqrt(a_L^3/mu);
n_orbits = 3; % Chosen display duration; Option A does not specify a duration.
tf = n_orbits*orb_period;

%% Orbital elements -> initial ECI states (Oct. 6 notes, p. 4)
[Lr0,Lv0] = oe2rv(a_L,e_L,i_L,bigW_L,wumg_L,nu_L0,mu);
[Fr0,Fv0] = oe2rv(a_F,e_F,i_F,bigW_F,wumg_F,nu_F0,mu);
leader0 = [Lr0;Lv0]; % 6-by-1 column, not a 3-by-2 array
follower0 = [Fr0;Fv0];

%% Initial relative state in the LEADER frame (Oct. 13 notes, p. 6)
C0 = eci2lvlh(Lr0,Lv0);
omega0 = [0;0;norm(cross(Lr0,Lv0))/norm(Lr0)^2];
rel_pos = C0*(Fr0-Lr0);
rel_vel = C0*(Fv0-Lv0)-cross(omega0,rel_pos);
x0 = [rel_pos;rel_vel];

fprintf('OPTION A: nonlinear two-body relative motion\n');
fprintf('Orbital period: %.6f s (%.6f min)\n',orb_period,orb_period/60);
fprintf('Simulation duration: %g periods\n',n_orbits);
fprintf('Initial rho [km]:       [% .9f % .9f % .9f]\n',rel_pos);
fprintf('Initial rho_dot [km/s]: [% .9f % .9f % .9f]\n',rel_vel);
fprintf('Initial separation: %.6f km\n',norm(rel_pos));

%% Simulink (the supplied, repaired OPTION A model)
load_system(fullfile(project_dir,'proj3_optiona.slx'));
simOut = sim('proj3_optiona','StopTime',num2str(tf,17), ...
    'ReturnWorkspaceOutputs','on');
out_rel = simOut.get('out_rel');
out_leader = simOut.get('out_leader');
time = out_rel.Time;
relative_state = out_rel.Data; % columns: x,y,z,xdot,ydot,zdot in LVLH
leader_state = out_leader.Data; % columns: r_ECI(3),v_ECI(3)
rho = relative_state(:,1:3);
rho_dot = relative_state(:,4:6);

%% Reconstruct follower ECI state from the simulated relative state
follower_state = zeros(size(leader_state));
for k = 1:numel(time)
    rL = leader_state(k,1:3).'; vL = leader_state(k,4:6).';
    C = eci2lvlh(rL,vL);
    omega = [0;0;norm(cross(rL,vL))/norm(rL)^2];
    follower_state(k,1:3) = (rL+C.'*rho(k,:).').';
    follower_state(k,4:6) = ...
        (vL+C.'*(rho_dot(k,:).'+cross(omega,rho(k,:).'))).';
end
separation = vecnorm(rho,2,2);
fprintf('Sampled separation range: %.6f to %.6f km\n', ...
    min(separation),max(separation));
fprintf('Final relative-position closure error: %.6g km\n', ...
    norm(rho(end,:).'-rel_pos));

%% Plots
labels = {'x (radial)','y (transverse)','z (normal)'};
figure('Name','Option A - Relative state','Color','w', ...
    'Position',[80 80 1080 720]);
for j = 1:3
    subplot(3,2,2*j-1);
    plot(time/60,rho(:,j),'LineWidth',1.4); grid on;
    ylabel([labels{j},' [km]']); xlabel('Time [min]');
    if j == 3, ylim([-1e-6,1e-6]); end % Do not magnify floating-point noise
    subplot(3,2,2*j);
    plot(time/60,rho_dot(:,j),'LineWidth',1.4); grid on;
    ylabel([labels{j},' velocity [km/s]']); xlabel('Time [min]');
    if j == 3, ylim([-1e-9,1e-9]); end
end
sgtitle('Follower relative to leader in rotating LVLH axes');
set(findall(gcf,'Type','axes'),'FontSize',11);

figure('Name','Option A - Relative path','Color','w');
plot(rho(:,1),rho(:,2),'LineWidth',1.5); hold on;
plot(rho(1,1),rho(1,2),'go','MarkerFaceColor','g');
plot(rho(end,1),rho(end,2),'rx','MarkerSize',9,'LineWidth',1.5);
axis equal; grid on;
xlabel('Radial x [km]'); ylabel('Transverse y [km]');
title('Relative path in the leader LVLH frame');
legend('Relative path','Start','End','Location','best');
set(gca,'FontSize',12);

figure('Name','Option A - ECI orbits','Color','w');
plot3(leader_state(:,1),leader_state(:,2),leader_state(:,3), ...
    'b-','LineWidth',1.5); hold on;
plot3(follower_state(:,1),follower_state(:,2),follower_state(:,3), ...
    'r--','LineWidth',1);
plot3(Lr0(1),Lr0(2),Lr0(3),'bo','MarkerFaceColor','b');
plot3(Fr0(1),Fr0(2),Fr0(3),'ro','MarkerFaceColor','r');
plot3(0,0,0,'k+','MarkerSize',10,'LineWidth',1.5);
axis equal; grid on; view(3);
xlabel('ECI X [km]'); ylabel('ECI Y [km]'); zlabel('ECI Z [km]');
title('Same orbit, different initial positions');
legend('Leader','Follower','Leader start','Follower start', ...
    'Earth center','Location','best');
set(gca,'FontSize',12);

figure('Name','Option A - Separation','Color','w');
plot(time/60,separation,'LineWidth',1.5); grid on;
xlabel('Time [min]'); ylabel('Separation [km]');
title('Leader-to-follower distance'); set(gca,'FontSize',12);

