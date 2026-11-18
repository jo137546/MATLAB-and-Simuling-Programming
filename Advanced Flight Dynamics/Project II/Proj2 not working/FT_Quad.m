function FT = FT_Quad(U)
phi = U(1); theta = U(2);

delta = [U(3);U(4);U(5);U(6)];

k1 = U(7); k2 = U(8); l = U(9);

mass = U(10); g = U(11);

FTM = [k1    k1      k1      k1;...
      0      0       0       l*k1;...
      l*k1 0      -l*k1  0;...
      -k2   k2      -k2     k2];
Prop_FT = FTM*delta;

F = Prop_FT(1);
tau_phi = Prop_FT(2);
tau_psi = Prop_FT(4);
torques = [tau_phi;tau_theta;tau_psi];


fb = [-g*sin(theta);
      g*cos(theta)*sin(phi);
      g*cos(theta)*cos(phi) - 1/mass*F];

FT = [fb;torques];