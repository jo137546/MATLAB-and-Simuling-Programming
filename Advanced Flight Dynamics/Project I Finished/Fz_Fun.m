function AOA = Fz_Fun(U)
m = U(1);
g = U(2);
rho = U(3);
V = U(4);
s = U(5);
beta = U(6);

fun = @(AOA)Fz_Balance(m,g,rho,V,s,AOA,beta);
AOA = fsolve(fun,5);
end

function F = Fz_Balance(m,g,rho,V,s,AOA,beta)
    F = 0.5*rho*V^2*s*((0.07+0.5*(-0.05+0.025*AOA)^2)*sind(AOA)*cosd(beta)+(-0.05+0.025*AOA)*cosd(AOA))-m*g;
end



