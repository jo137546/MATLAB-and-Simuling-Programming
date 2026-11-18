%{
EML 6405 Project 1
Elliott Eggers
2/15/23
%}
clear;clc;close all

%% Constant Parameters
m = 2; % kg
g = 9.81; % m/s^2
p = 0;
q = 0;
r = 0;
theta = 0;
phi = 0;
psi = 0;
Re = 6.3781*10^6; % m (Radius of Earth)

Isp = 50;
Ca = 5;
Ac = pi*0.125^2;

beta = 0;
s = 0.152*0.762;

t_final = 100;
t_step = 0.01;
time = 0:t_step:t_final;

%% Initial Conditions
alt_0 = 5000; %try to match xu notes

pos_0 = [0;0;-alt_0-Re];
rho_0 = density(-pos_0(3)-Re);

Vel_0 = [30;0;0];
V_0 = sqrt(Vel_0(1)^2+Vel_0(2)^2+Vel_0(3)^2);

%U = [m,g,rho_0,V_0,s,beta];
%AOA = Fz_Fun(U);

[AOA,del_t,theta,Vel_0] = SLUF_Trim(m,g,rho_0,V_0,s,Isp,Ca,Ac); %for SLUF Flight

C = lift_drag_fun(AOA);
CL = C(1);
CD = C(2);

FL = 0.5*rho_0*V_0^2*s*-CL;
FD = 0.5*rho_0*V_0^2*s*-CD;

F = T_W2B(AOA,beta)*[FD;0;FL];

%del_t = -F(1)/(0.029*Isp*V_0*rho_0*g*Ac*Ca);

%% Simulation
sim("Project1_Sim.slx")

%% Plot Results
X = ans.X;
Y = ans.Y;
Z = -ans.Z-Re;

figure('color','w')
subplot(3,1,1)
plot(time,X,'k')
xlabel('time [s]')
ylabel('North [m]')

subplot(3,1,2)
plot(time,Y,'k')
xlabel('time [s]')
ylabel('East [m]')

subplot(3,1,3)
plot(time,Z,'k')
xlabel('time [s]')
ylabel('Altitude [m]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    A/C visualizer     &&&&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trajectory and compass heading
heading = mod(psi,360);

figure('Color','w')

% Flight path: horizontal axes East/North, vertical axis altitude
subplot(1,2,1)
plot3(Y,X,Z,'b','LineWidth',1.5)
hold on
plot3(Y(end),X(end),Z(end),'ro','MarkerFaceColor','r')
grid on
xlabel('East [m]')
ylabel('North [m]')
zlabel('Altitude [m]')
title('SLUF trajectory')
view(3)

% Compass: North is up, East is right
subplot(1,2,2)
hold on
axis equal
axis([-1.3 1.3 -1.3 1.3])
axis off

a = linspace(0,360,361);
plot(sind(a),cosd(a),'k','LineWidth',1.2)

text(0, 1.15,'N','HorizontalAlignment','center')
text(1.15,0,'E','HorizontalAlignment','center')
text(0,-1.15,'S','HorizontalAlignment','center')
text(-1.15,0,'W','HorizontalAlignment','center')

% Arrow components: East = sin(heading), North = cos(heading)
quiver(0,0,0.85*sind(heading),0.85*cosd(heading),0, ...
    'r','LineWidth',2,'MaxHeadSize',0.4)

title(sprintf('Heading: %.1f deg',heading))

%% =========================================================
%  COORDINATE FRAME / ATTITUDE VISUALIZATION
% =========================================================

% SLUF aircraft attitude
phi   = 0;       % roll angle [deg]
theta = AOA;     % pitch angle [deg]
psi   = 0;       % yaw angle [deg]

% ---------------------------------------------------------
% Principal coordinate transformations
% ---------------------------------------------------------

C_x = CoordinateTransform('x',phi);
C_y = CoordinateTransform('y',theta);
C_z = CoordinateTransform('z',psi);

% NED -> Body transformation
C_N2B = C_x * C_y * C_z;

% Body -> NED transformation
% Rotation matrices are orthogonal, so inverse = transpose
C_B2N = C_N2B.';

% ---------------------------------------------------------
% Aircraft body axes expressed in NED coordinates
% ---------------------------------------------------------

xB = C_B2N * [1;0;0];
yB = C_B2N * [0;1;0];
zB = C_B2N * [0;0;1];

% ---------------------------------------------------------
% Plot coordinate systems
% ---------------------------------------------------------

figure('Color','w')
hold on
grid on
axis equal

L = 1;

% NED frame
quiver3(0,0,0,L,0,0,'LineWidth',2)
quiver3(0,0,0,0,L,0,'LineWidth',2)
quiver3(0,0,0,0,0,L,'LineWidth',2)

text(L,0,0,'N')
text(0,L,0,'E')
text(0,0,L,'D')

% Aircraft body frame
quiver3(0,0,0,...
    L*xB(1),L*xB(2),L*xB(3),...
    'LineWidth',2)

quiver3(0,0,0,...
    L*yB(1),L*yB(2),L*yB(3),...
    'LineWidth',2)

quiver3(0,0,0,...
    L*zB(1),L*zB(2),L*zB(3),...
    'LineWidth',2)

text(L*xB(1),L*xB(2),L*xB(3),'x_B')
text(L*yB(1),L*yB(2),L*yB(3),'y_B')
text(L*zB(1),L*zB(2),L*zB(3),'z_B')

xlabel('North')
ylabel('East')
zlabel('Down')

set(gca, 'ZDir', 'reverse')

title(sprintf(...
    'SLUF Attitude: \\alpha = %.2f deg, \\theta = %.2f deg',...
    AOA,theta))

view(135,25)
axis equal
grid on