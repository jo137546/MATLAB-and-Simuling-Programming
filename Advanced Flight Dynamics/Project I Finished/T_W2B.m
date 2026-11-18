%{
Readme
Written by: Elliott Eggers
1/25/2023

Creates Transformation Matrix from Wind Coordinates to Body Coordinates by
angle of attack (a) and sideslip angle (B)
%}

function R = T_W2B(a,B)
R = [cosd(a) 0 -sind(a); 0 1 0; sind(a) 0 cosd(a)]*[cosd(B) -sind(B) 0; sind(B) cosd(B) 0; 0 0 1];
end