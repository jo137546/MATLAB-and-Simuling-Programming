%{
Readme
Written by Elliott Eggers
1/30/23

Approximates the coefficients of lift and drag as a function of the angle
of attack (a) in degrees
%}


function C = lift_drag_fun(a)
CL = -0.05+0.025*a;
CD = 0.07+0.5*CL^2;

C = [CL,CD];
end