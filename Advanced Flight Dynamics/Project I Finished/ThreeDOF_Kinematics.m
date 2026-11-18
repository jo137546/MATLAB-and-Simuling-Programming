%{
Read me
Written by Elliott Eggers 
1/30/23

Calculates dr/dt in NED coordinates from [u,v,w]' in body coordinates using
the Euler angles in degrees. drdt is a vector which contains [dx/dt dy/dt
dz/dt]. The input vector U contains [u,v,w,theta,phi,psi)]
%}

function drdt = ThreeDOF_Kinematics(U)
u = U(1);
v = U(2);
w = U(3);
theta = U(4);
phi = U(5);
psi = U(6);

%corrected terms
drdt = [cosd(theta)*cosd(psi) sind(phi)*sind(theta)*cosd(psi)-cosd(phi)*sind(psi) cosd(phi)*sind(theta)*cosd(psi)+sind(phi)*sind(psi);...
    cosd(theta)*sind(psi) sind(phi)*sind(theta)*sind(psi)+cosd(phi)*cosd(psi) cosd(phi)*sind(theta)*sind(psi)-sind(phi)*cosd(psi);...
    -sind(theta) sind(phi)*cosd(theta) cosd(phi)*cosd(theta)]*[u;v;w];
end