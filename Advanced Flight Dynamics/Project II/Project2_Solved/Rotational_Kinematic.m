function y = Rotational_Kinematic(u)

pqr = [u(1);u(2);u(3)]; %[p,q,r] roll, pitch, yaw angles
phi = u(4);
theta = u(5);

Rot_Kine_M =  [1    sin(phi)*tan(theta) cos(phi)*tan(theta);...
               0    cos(phi)          -sin(phi);...
               0    sin(phi)/cos(theta) cos(phi)/cos(theta)];
           
y = Rot_Kine_M*pqr;  %[phi_dot,theta_dot,psi_dot]