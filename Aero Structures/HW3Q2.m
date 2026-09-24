n = linspace(2, 200, 1000);
t = 1;
a = n*t;

Ihollow = ( 2 * (a.^3)) ./ 3;
Isolid = ( (a.^ 4) ) ./ 12;
Iratio = Ihollow./Isolid;

plot(a , Iratio, 'r')

xlabel('a')
ylabel('Moment of Inertia Ratio')
title('Plot of Moment of inertia ratio')
legend('Iratio')

hold off