
function y = Trans_Kinematic(U)

uvw = [U(1);U(2);U(3)];   %u, v ,w

phi = U(4); theta = U(5); psi = U(6);

Trans_Kine_M = [cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);...
                cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);...
                sin(theta)    -sin(phi)*cos(theta) -cos(phi)*cos(theta)];


