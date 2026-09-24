% Laplace / Inverse Laplace Interactive Tool
%clear; clc;
syms t s

% Ask user for choice
choice = input('Enter 1 for Laplace Transform, 2 for Inverse Laplace: ');

if choice == 1
    % Laplace Transform
    disp('Enter a time-domain function f(t). Use variable t.');
    f = input('f(t) = ');
    
    F = laplace(f, t, s);
    
    disp('Laplace Transform F(s) = ');
    pretty(F)
    
elseif choice == 2
    % Inverse Laplace Transform
    disp('Enter an s-domain function F(s). Use variable s.');
    F = input('F(s) = ');
    
    f = ilaplace(F, s, t);
    
    disp('Inverse Laplace f(t) = ');
    pretty(f)
    
else
    disp('Invalid choice. Please enter 1 or 2.');
end
