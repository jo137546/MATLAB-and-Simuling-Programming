function qs = qs_Fun(U)
m = U(1);
g = U(2);
AOA = U(3);
beta = U(4);
CL = U(5);
CD = U(6);



%qs = (m*g)/(CD*sind(AOA)*cosd(beta)+CL*cosd(AOA)); incorrect!

qs = (m*g*cosd(AOA)) / ...
     (CD*sind(AOA)*cosd(beta) + CL*cosd(AOA)); % COMMAND PRESSURE not dynamic pressure model 
                                               % becuase theta = AOA body weight is mg cos(AOA)

