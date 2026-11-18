%{
Read me
Written by Elliott Eggers
1/30/23

Calculates the thrust force applied to the A/C from a turbojet, ramjet,
scramjet, or turbofan in Mach <3 U = [delta,Isp,g,rho,M,a,Ca,Ac]
%}

function Ft = thrust(U)
delta = U(1);
Isp = U(2);
g = U(3);
rho = U(4);
M = U(5);
a = U(6);
Ca = U(7);
Ac = U(8);

Ft = 0.029*delta*Isp*g*rho*M*a*Ca*Ac;
end