%Declare any needed variables
syms M_0 L b EI z a o C_1 C_2

% Define the support reactions shown in Fig. 16.23
R_A=-M_0/L;
R_C=M_0/L;

% Define the equation for the bending moment at any section Z between B and C
M=[-R_A*z -M_0*(a)^o]; % For z<¼b: a¼0; For z>b: a¼z-b

% Substitute M into the second of Eq. (16.32)
v_zz=-M/EI;

% Integrate v_zz to get v_z and v
v_z=[int(v_zz(1),z)+C_1/EI int(v_zz(2),a)];
v=[int(v_z(1),z)+C_2/EI int(v_z(2),a)];
v=sum(subs(v,o,0));

% Use boundary conditions to determine C_1 and C_2
% BC #1: v¼0 when z¼0
c_2=solve(subs(subs(v*EI,a,0),z,0),C_2);
v=subs(v,C_2,c_2);

% BC #1: v¼0 when z¼0
c_1=solve(subs(subs(v*EI,a,z-b),z,L),C_1);
v=simplify(subs(v,C_1,c_1));

% Output the resulting deflection equation to the Command Window
disp('The equation for the deflection curve of the beam is:')
disp(['v=' char(v)])
disp('Where: a=0 for z<=b, and a=z-b for z>b')