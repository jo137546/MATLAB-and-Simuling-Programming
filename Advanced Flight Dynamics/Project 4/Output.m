function Y = Output(u)
x = u(1:6);
[~,~,C] = CW_Matrix;
Y = C*x;

%Input = [6;1]
%Output = [3;1]


