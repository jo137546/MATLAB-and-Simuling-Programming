clear; clc; close;

%%
global J k1 k2 l mt g delta x0 a0 vel_initial
k1 = 1; k2 = 0.75;
M = 0.25; m = 0.01; mt = M+4*m;
l = 0.2; R = 0.03;
J = MomentsofInertia([M,R,m,l]);
g = 9.8;

delta_hov = (mt*g/4*k1)*ones(4,1);

T_phi0 = 0.05;
drl_yaw = (mt*g/4*k1) + (T_phi0/4*k2);
delta_yaw = [dfb_yaw;drl_yaw;dfb_yaw;drl_yaw];

delta = delta_yaw;
%%
x0 = [0;0;10]; a0 = [0;0;0];

v0 = [0;0;0]; da0 = [0;0;0];

vel_initial = [v0;da0];
%%
tf = 15;
%%
sim('Project2_Sim_V2.slx');

%%
figure(1);
subplot(1,2,1);
plot(ans.POSITION,'LineWidth',1.5);
axis([0 15 -5 15])
xlabel('Time (s)','fontweight','bold'); ylabel('Position (m)','fontweight','bold');
legend('x','y','z');
title('Position');

subplot(1,2,2);
plot(ans.PHI,'LineWidth',1.5); hold on; plot(ans.THETA,'LineWidth',1.5);
hold on; plot(ans.PSI,'LineWidth',1.5);
axis([0 15 -5 5])
legend('\phi','\theta','\psi'); title('Direction Angle');
xlabel('Time (s)','fontweight','bold'); ylabel('Direction Angle (deg)','fontweight','bold');
