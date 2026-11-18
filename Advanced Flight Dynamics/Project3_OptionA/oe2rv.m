function [r_bar,v_bar] = oe2rv(a,e,i,bigW,wumg,nu,mu)
% Classical orbital elements -> ECI position and velocity.
% a in km; i, RAAN (bigW), argument of periapsis (wumg), nu in radians.
% Oct. 6 spacecraft notes, p. 4; equivalent perifocal rotation form.
if nargin < 7, mu = 3.986e5; end
p = a*(1-e^2);
r_pf = p/(1+e*cos(nu))*[cos(nu);sin(nu);0];
v_pf = sqrt(mu/p)*[-sin(nu);e+cos(nu);0];

R3_Omega = [cos(bigW),-sin(bigW),0;sin(bigW),cos(bigW),0;0,0,1];
R1_i = [1,0,0;0,cos(i),-sin(i);0,sin(i),cos(i)];
R3_omega = [cos(wumg),-sin(wumg),0;sin(wumg),cos(wumg),0;0,0,1];
Q = R3_Omega*R1_i*R3_omega; % Perifocal -> ECI
r_bar = Q*r_pf;
v_bar = Q*v_pf;
end

