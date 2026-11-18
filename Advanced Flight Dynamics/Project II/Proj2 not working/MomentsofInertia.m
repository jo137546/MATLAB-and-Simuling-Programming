function J = MomentsofInertia(U)
M = U(1); R = U(2);
m = U(3); l = U(4);

Jx = (2/5)*M*R^2+2*l^2*m;
Jy = Jx;
Jz = (2/5)*M*R^2+4*l^2*m;
J = [Jx;Jy;Jz];