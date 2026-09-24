clear; clc;

% --- matrices (numeric doubles) ---
A = [2 -2 -5/2;
     0 -8 -2;
     0 0 -7];

B = [1; 1; 1;];
C = [1 0 1];
D = 0;

A = double(A); B = double(B); C = double(C); D = double(D);

% Build SS and TF
sys_ss = ss(A,B,C,D);
sys_tf = tf(sys_ss);

disp('Transfer function G(s) = Y(s)/U(s):')
sys_tf

% Given sys_ss = ss(A,B,C,D);
G = tf(sys_ss);

% numerical noise and scale to integers
[num,den] = tfdata(G,'v');
tol = 1e-12;
num(abs(num)<tol) = 0;
den(abs(den)<tol) = 0;

% scale to match the earlier integer form
scale = 1;                     % choose 2 to get integer coefficients
num = scale*num; den = scale*den;

G_clean = tf(num,den)
