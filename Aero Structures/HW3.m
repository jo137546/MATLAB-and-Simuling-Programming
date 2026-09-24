n = linspace(2, 200, 1000);
t = 0.1;
a = n*t;

Iapprox = ( 2 * (a.^3) * t) / 3;
Iactual = ( (a.^4) - (a - (2 * t) ).^4) / 12;

plot(a,Iapprox, 'r')
hold on
plot(a,Iactual, 'b')

error = abs((Iactual - Iapprox) ./ Iactual) * 100;
max_error = max(error);
min_error = min(error);

xlabel('a')
ylabel('Moment of Inertia')
title('Plot of the Sine Function')
legend('Iapprox', 'Iactual')

disp(['Maximum error: ', num2str(max_error), '%'])
disp(['Minimum error: ', num2str(min_error), '%'])