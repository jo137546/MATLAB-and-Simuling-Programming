function y = leaderP(u,mu)
% Pack one output for the Simulink Interpreted MATLAB Function block.
% y(1:3): leader's Earth-centered position expressed in its LVLH axes.
% y(4): leader true anomaly (rad).
if nargin < 2, mu = 3.986e5; end
[~,~,~,~,~,nu] = orbit_ele(u,mu);
rL_vec = [norm(u(1:3));0;0]; % Oct. 13 relative-motion notes, p. 3
y = [rL_vec;nu];
end

