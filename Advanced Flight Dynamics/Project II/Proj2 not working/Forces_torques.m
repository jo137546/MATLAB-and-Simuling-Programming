function y = Forces_torques(u)
deltas = [u(1);u(2);u(3);u(4)]; %[delta_f;delta_r;delta_b;delta_l]
phi = u(5);
theta = u(6);
kt = u(7);  %thrust coefficient k1
kd = u(8);  %drag coefficient k2
arm = u(9);
mm = u(10);
g = 9.8;

MM = [kt    kt      kt      kt;...
      0     0       0       arm*kt;...
      arm*kt 0      -arm*kt  0;...
      -kd   kd      -kd     kd];
temp = MM*deltas;
F = temp(1);
tau_phi = temp(2);
tau_theta = temp(3);
tau_psi = temp(4);
torques = [tau_phi;tau_theta;tau_psi];

fx = -g*sin(theta);
fz = g*cos(theta)*cos(phi) - 1/mm*F;
fb = [fx;fy;fz];

y = [fb;torques];