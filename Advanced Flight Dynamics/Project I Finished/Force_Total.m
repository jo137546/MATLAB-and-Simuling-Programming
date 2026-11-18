%{
Readme:
Written by Elliott Eggers
2/6/23

Calculates the total force vector from the aerodynamic forces and the
thrust U = [CL,CD,qs,alpha,beta,Ft] and F = [Fx,Fy,Fz]
%}

function F = Force_Total(U)

CL = U(1);
CD = U(2);
qs = U(3);
alpha = U(4);
beta = U(5);
Ft = U(6);


Fa = [cosd(alpha) 0 -sind(alpha);0 1 0;sind(alpha) 0 cosd(alpha)]*...
     [cosd(beta) -sind(beta) 0;sind(beta) cosd(beta) 0;0 0 1]*...
     [-qs*CD; 0; -qs*CL];

F(1) = Fa(1)+Ft;
F(2) = Fa(2);
F(3) = Fa(3);

end