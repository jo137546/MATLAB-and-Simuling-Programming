function y = leader_angV(u,a,e,mu)
% Leader LVLH angular rate and angular acceleration about local z.
% Oct. 13 relative-motion notes, p. 3 (elliptic-orbit case).
if nargin < 4, mu = 3.986e5; end
if nargin < 3, e = 0.05; end
if nargin < 2, a = 8000; end
nu = u(1);
p = a*(1-e^2);
r = p/(1+e*cos(nu));
h = sqrt(mu*p);
wumg_nu = h/r^2; % nu_dot, rad/s
r_dot = mu/h*e*sin(nu);
acc_nu = -2*wumg_nu*r_dot/r; % nu_ddot, rad/s^2
y = [wumg_nu;acc_nu];
end

