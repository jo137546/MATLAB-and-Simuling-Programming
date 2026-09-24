n = linspace(2, 200, 1000);
t = 0.1;
a = n*t;

Iapprox = ( 2 * (a.^3) * t) / 3;
Iactual = ( (a.^4) - (a - (2 * t) ).^4) / 12;

plot(a,Iapprox, 'r')
hold on
plot(a,Iactual, 'b')

error = abs((Iactual - Iapprox) ./ Iactual) * 100;
[~, idx_10] = min(abs(a/t - 10)); % Find the closest value of a/t = 10
[~, idx_25] = min(abs(a/t - 25)); % Find the closest value of a/t = 25
[~, idx_100] = min(abs(a/t - 100)); % Find the closest value of a/t = 100

error_10 = error(idx_10);   % Error at a/t = 10
error_25 = error(idx_25);   % Error at a/t = 25
error_100 = error(idx_100); % Error at a/t = 100

xlabel('a')
ylabel('Moment of Inertia')
title('Plot of the Sine Function')
legend('Iapprox', 'Iactual')

disp(['Error at a/t = 10: ', num2str(error_10), '%'])
disp(['Error at a/t = 25: ', num2str(error_25), '%'])
disp(['Error at a/t = 100: ', num2str(error_100), '%'])