function dy = dynamics(U)

u = U(1); v = U(2); w = U(3);

p = U(4); q = U(5); r = U(6);

fb = [U(7);U(8);U(9)];
T_phi = U(10); T_theta = U(11); T_psi = U(12);

Jx = U(13); Jy = U(14); Jz = U(15);
%Translational dynamics
d_uvw = fb + [r*v-q*w; p*w-r*u; q*u-p*v];

d_pqr = [(T_phi + (Jy-Jz)*q*r)/Jx;
         (T_theta + (Jz-Jx)*p*r)/Jy;
         (T_psi + (Jx-Jy)*p*q)/Jz];

dy = [d_uvw;d_pqr];
