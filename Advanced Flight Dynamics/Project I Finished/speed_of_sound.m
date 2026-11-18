%{
Readme:
Written by Elliott Eggers
1/14/2023

Calculates the speed of sound in air given the temperature.

In metric, speed is given in meters/second

%}


function a = speed_of_sound(T)
R = 287.052874; %J/(kg*K)
T = T+273.15;

gamma = 1.4;

a = sqrt(gamma*R*T);
end
