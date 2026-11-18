function dy = rel_dyn(u,mu)
% Exact nonlinear relative dynamics in the leader's rotating LVLH frame.
% Oct. 13 relative-motion notes, pp. 2-3. No control/perturbation forces.
% Mux: [nu_L; nu_dot_L; nu_ddot_L; rho(3); rho_dot(3); rL_LVLH(3)].
if nargin < 2, mu = 3.986e5; end
leader_wumg = [0;0;u(2)];
leader_acc = [0;0;u(3)];
rol = u(4:6); rol = rol(:);
rol_v = u(7:9); rol_v = rol_v(:);
rL_vec = u(10:12); rL_vec = rL_vec(:);
rF_vec = rL_vec+rol;
rL_mag = norm(rL_vec);
rF_mag = norm(rF_vec);

rol_acc = -mu*rF_vec/rF_mag^3 + mu*rL_vec/rL_mag^3 ...
    -2*cross(leader_wumg,rol_v) ...              % Coriolis
    -cross(leader_acc,rol) ...                   % Euler acceleration
    -cross(leader_wumg,cross(leader_wumg,rol));  % Centrifugal
dy = [rol_v;rol_acc];
end

