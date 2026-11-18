function dx = CW_Equation(u)
x = u(1:6);
U = u(7:9);

[A,B] = CW_Matrix;

dx = A*x + B*U;

%Input = [9;1]
%Output = [6;1]