
% Clear command window, close all graphs, clear workspace
clc, close all, clear

% Define the function f(x)
f = @(x) -0.5-2.25*(x*cosh(1.5*x)) ;      % Insert Function f(x)
df = @(x) -2.25*cosh(1.5*x)-2.25*x*((1.5)*(sinh(1.5*x)));     % Insert Exact First Derivative df(x)
% Select x value to evaluate derivative
x = 1.75;
% loop through the exponents for each delta x value and evaluate the derivative
for i = 1:20
  
  % calculate the delta x
  del(i,1) = 10^(-i);
  % Backward Finite Difference
  backward(i,1) = (f(x)-f(x-del(i)))/del(i);
  % Forward Finite Difference
  forward(i,1) = (-f(x)+f(x+del(i)))/del(i);
  % Central Finite Difference
  central(i,1) = (f(x+del(i))-f(x-del(i)))/2*del(i);
  % calculate the errors compared to exact (analytical) derivative
  backward_error(i,1) = abs(backward(i)-df(x));
  forward_error(i,1) = abs(forward(i)-df(x));
  central_error(i,1) = abs(central(i)-df(x));
  
end
% Plot the error vs del x for each method on the log scale
figure
hold on
% the loglog() function is the same as plot() but on log scale for both axes
loglog(del,backward_error)
loglog(del,forward_error)
loglog(del,central_error)
xlabel('\Delta'), ylabel('Error')
title('Error on Fordward, Backward, and Central')
legend('Backward','Forward','Central')
grid on, hold off
% assemble the results into a table (actually a matrix)
Table = [del,backward_error,forward_error,central_error];


% print to screen
print("=========== Double Precision Results ===========\n")
disp(Table)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repeat everything, but in single precision
% Select x value to evaluate derivative
x = single(x);
for i = single(1:20)
  
end
% Plot the error vs del x for each method on the log scale
figure
hold on
% the loglog() function is the same as plot() but on log scale for both axes
loglog(del,backward_error)
loglog(del,forward_error)
loglog(del,central_error)
xlabel('\Delta'), ylabel('Error')
title('Error on Fordward, Backward, and Central')
legend('Forward','Backward','Central')
grid on
% assemble the results into a table (actually a matrix)
fprintf("\n\n\n")
Table = [];
% print to screen
printf("=========== Single Precision Results ===========\n")
disp(Table);