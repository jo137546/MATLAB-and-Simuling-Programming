function dy = lead_orb(u,mu)
% Two-body leader dynamics in ECI. Units: km, km/s, s.
% Oct. 1 spacecraft notes, pp. 1 and 4.
if nargin < 2, mu = 3.986e5; end
r_bar = u(1:3);
v_bar = u(4:6);
dy = [v_bar(:); -mu*r_bar(:)/norm(r_bar)^3];
end

