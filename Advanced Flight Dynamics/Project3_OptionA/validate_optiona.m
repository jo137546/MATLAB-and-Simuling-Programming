% Run AFTER proj3a_main. Checks the nonlinear Simulink output against an
% independent analytic two-body solution obtained from Kepler's equation.
assert(exist('relative_state','var') == 1,'Run proj3a_main first.');
assert(all(isfinite(relative_state),'all'),'Non-finite simulation output.');
assert(numel(time)==size(leader_state,1),'State time grids do not match.');
assert(max(abs(time-out_leader.Time)) < 1e-10,'State time grids do not match.');

% Orbital-element round trip: includes the negative initial true anomaly.
for initial_state = [leader0,follower0]
    [aa,ee,ii,OO,ww,ff] = orbit_ele(initial_state,mu);
    [rr,vv] = oe2rv(aa,ee,ii,OO,ww,ff,mu);
    assert(norm(rr-initial_state(1:3)) < 1e-8,'Position round trip failed.');
    assert(norm(vv-initial_state(4:6)) < 1e-11,'Velocity round trip failed.');
end
assert(norm(C0*C0.'-eye(3),'fro') < 1e-12 && det(C0)>0, ...
    'LVLH basis is not a proper orthonormal frame.');

% Closed-form reference: advance mean anomaly, solve Kepler's equation,
% then convert the eccentric anomaly directly to Cartesian position.
leader_exact = kepler_reference(time,a_L,e_L,i_L,bigW_L,wumg_L,nu_L0,mu);
follower_exact = kepler_reference(time,a_F,e_F,i_F,bigW_F,wumg_F,nu_F0,mu);
relative_exact = zeros(size(relative_state));
for k = 1:numel(time)
    r = leader_exact(k,1:3).'; v = leader_exact(k,4:6).';
    C = eci2lvlh(r,v);
    om = [0;0;norm(cross(r,v))/norm(r)^2];
    rp = C*(follower_exact(k,1:3).'-r);
    rv = C*(follower_exact(k,4:6).'-v)-cross(om,rp);
    relative_exact(k,:) = [rp;rv].';
end
position_error = vecnorm(relative_state(:,1:3)-relative_exact(:,1:3),2,2);
velocity_error = vecnorm(relative_state(:,4:6)-relative_exact(:,4:6),2,2);
max_pos_error = max(position_error);
max_vel_error = max(velocity_error);
max_z = max(abs(rho(:,3)));
closure_error = norm(rho(end,:).'-rel_pos);

fprintf('\nVALIDATION AGAINST ANALYTIC KEPLER SOLUTION\n');
fprintf('Max relative-position error: %.9g km\n',max_pos_error);
fprintf('Max relative-velocity error: %.9g km/s\n',max_vel_error);
fprintf('Max out-of-plane displacement: %.9g km\n',max_z);
fprintf('Relative-position closure after %g periods: %.9g km\n', ...
    n_orbits,closure_error);
assert(max_pos_error < 1e-3,'Relative-position error exceeds 1 m.');
assert(max_vel_error < 1e-6,'Relative-velocity error exceeds 1 mm/s.');
assert(max_z < 1e-6,'Unexpected out-of-plane motion.');
assert(closure_error < 1e-3,'Relative orbit did not close within 1 m.');

states = {leader_state,follower_state};
for j = 1:2
    state = states{j};
    rmag = vecnorm(state(:,1:3),2,2);
    energy = sum(state(:,4:6).^2,2)/2-mu./rmag;
    hvec = cross(state(:,1:3),state(:,4:6),2);
    energy_drift = max(abs(energy-energy(1)))/abs(energy(1));
    h_drift = max(vecnorm(hvec-hvec(1,:),2,2))/norm(hvec(1,:));
    fprintf('Satellite %d: fractional energy drift %.9g, angular-momentum drift %.9g\n', ...
        j,energy_drift,h_drift);
    assert(energy_drift < 1e-8 && h_drift < 1e-8, ...
        'Two-body invariant conservation check failed.');
end
fprintf('PASS: all Option A validation checks.\n');

function state = kepler_reference(t,a,e,i,Omega,w,nu0,mu)
n = sqrt(mu/a^3);
E0 = atan2(sqrt(1-e^2)*sin(nu0),e+cos(nu0));
M = E0-e*sin(E0)+n*t(:);
E = M;
for j = 1:12
    E = E-(E-e*sin(E)-M)./(1-e*cos(E));
end
assert(max(abs(E-e*sin(E)-M)) < 1e-11,'Kepler solver failed.');
Q = [cos(Omega)*cos(w)-sin(Omega)*sin(w)*cos(i), ...
    -cos(Omega)*sin(w)-sin(Omega)*cos(w)*cos(i),sin(Omega)*sin(i); ...
    sin(Omega)*cos(w)+cos(Omega)*sin(w)*cos(i), ...
    -sin(Omega)*sin(w)+cos(Omega)*cos(w)*cos(i),-cos(Omega)*sin(i); ...
    sin(w)*sin(i),cos(w)*sin(i),cos(i)];
rpf = [a*(cos(E)-e),a*sqrt(1-e^2)*sin(E),zeros(size(E))];
vpf = [-a*n*sin(E)./(1-e*cos(E)), ...
    a*n*sqrt(1-e^2)*cos(E)./(1-e*cos(E)),zeros(size(E))];
state = [rpf*Q.',vpf*Q.'];
end

