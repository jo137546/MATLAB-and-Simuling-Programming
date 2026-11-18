function y = Rotational_kine(u)

RPY = [u(1);u(2);u(3)]; %[p,q,r] roll, pitch, yaw angles
phi = u(4);
R = RMatrixT_pqr2rates(phi,theta);
y = R'*RPY;  %[phi_dot,theta_dot,psi_dot]