%{
Readme
Written by Elliott Eggers
2/4/23

Calculates 3DOF Kinetics for an A/C in body coordinates. Here rdot is a
vector which contains udot, vdot, and wdot. The input U is a vector which
contains the following variables: [m,g,p,q,r,u,v,w,theta,phi,Fx,Fy,Fz]
%}

function rdot = ThreeDOF_Kinetics(U)
m = U(1);
g = U(2);
p = U(3);
q = U(4);
r = U(5);
u = U(6);
v = U(7);
w = U(8);
theta = U(9);
phi = U(10);
Fx = U(11);
Fy = U(12);
Fz = U(13);

rdot = [Fx/m-g*sind(theta)-q*w+r*v;...
        Fy/m+g*cosd(theta)*sind(phi)-r*u+p*w;...
        Fz/m+g*cosd(theta)*cosd(phi)-p*v+q*u];
end