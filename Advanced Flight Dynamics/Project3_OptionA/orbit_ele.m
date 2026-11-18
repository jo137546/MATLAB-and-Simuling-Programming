function [a,e,i,bigW,wumg,nu] = orbit_ele(u,mu)
% ECI state -> classical orbital elements. Units: km and radians.
% Oct. 6 spacecraft notes, pp. 5-6. atan2 keeps angle quadrants correct.
% This project uses a noncircular, nonequatorial elliptic orbit.
if nargin < 2, mu = 3.986e5; end
r_bar = u(1:3); r_bar = r_bar(:);
v_bar = u(4:6); v_bar = v_bar(:);
r = norm(r_bar);
h = cross(r_bar,v_bar); H = norm(h);
n = cross([0;0;1],h); N = norm(n);
e_bar = cross(v_bar,h)/mu-r_bar/r;
e = norm(e_bar);
a = 1/(2/r-dot(v_bar,v_bar)/mu);
i = atan2(norm(h(1:2)),h(3));
if e < 1e-10 || N/H < 1e-10
    error('orbit_ele:DegenerateOrbit', ...
        'This helper requires a noncircular, nonequatorial orbit.');
end
bigW = mod(atan2(n(2),n(1)),2*pi);
wumg = mod(atan2(dot(cross(n,e_bar),h)/(N*e*H), ...
    dot(n,e_bar)/(N*e)),2*pi);
nu = atan2(dot(cross(e_bar,r_bar),h)/(e*r*H), ...
    dot(e_bar,r_bar)/(e*r)); % Signed: [-pi,pi], including initial -10 deg
end

