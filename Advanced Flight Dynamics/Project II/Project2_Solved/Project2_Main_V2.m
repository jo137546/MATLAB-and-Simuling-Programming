clear; clc; close;

%%
global J k1 k2 l mt g delta x0 a0 vel_initial
k1 = 1; k2 = 0.75;
M = 0.25; m = 0.01; mt = M+4*m;
l = 0.2; R = 0.03;
J = MomentsofInertia([M,R,m,l]);
g = 9.8;

delta_hov = (mt*g/(4*k1))*ones(4,1);
delta = delta_hov;
tf = 15;
flight_cases = 1:3; % Use 1, 2, or 3 here to run only one case.
case_names = {'Alt Hold','Alt Hold + Spin','SLUF'};

for flight_case = flight_cases
%%
x0 = [0;0;10]; a0 = [0;0;0];

v0 = [0;0;0]; da0 = [0;0;0];
if flight_case == 2
    da0(3) = 20*pi/180; % Constant yaw rate, rad/s.
elseif flight_case == 3
    v0(1) = 5; % Straight level flight, m/s.
end

vel_initial = [v0;da0];
%%
simOut = sim('Project2_Sim_V2.slx');

%%
figure(flight_case); clf;
subplot(1,2,1);
plot(simOut.POSITION,'LineWidth',1.5);
xlim([0 tf])
xlabel('Time (s)','fontweight','bold'); ylabel('Position (m)','fontweight','bold');
legend('x_N','y_E','h');
title([case_names{flight_case} ': Position']);

subplot(1,2,2);
plot(simOut.PHI.Time,simOut.PHI.Data*180/pi,'LineWidth',1.5); hold on;
plot(simOut.THETA.Time,simOut.THETA.Data*180/pi,'LineWidth',1.5);
plot(simOut.PSI.Time,mod(simOut.PSI.Data,2*pi)*180/pi,'LineWidth',1.5);
xlim([0 tf])
legend('\phi','\theta','\psi'); title('Direction Angle');
xlabel('Time (s)','fontweight','bold'); ylabel('Direction Angle (deg)','fontweight','bold');

saveas(gcf,sprintf('Project2_Case%d.png',flight_case));
save(sprintf('Project2_Case%d.mat',flight_case),'simOut','delta','x0','a0','v0','da0','tf');
end
