function dy = dynamics(U)

u = U(1); v = U(2); w = U(3);

p = U(4); q = U(5); r = U(6);

fb = [U(7);U(8);U(9)];
T_phi = U(10); T_theta = U(11); T_psi = U(12);

Jx = U(13); Jy = U(14); Jz = U(15);
%Translational dynamics
d_uvw = fb;

dy = [d_uvw;d_pqr];
