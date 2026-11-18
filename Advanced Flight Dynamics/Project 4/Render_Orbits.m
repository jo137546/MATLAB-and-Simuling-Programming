function htmlFile = Render_Orbits(time, state, target)
% Render_Orbits  Export the simulated motion to an offline 3-D spacecraft scene.
% Usage after a simulation: Render_Orbits(time, x, rho_f)
% time: N-by-1 seconds; state: N-by-6 [x y z vx vy vz], LVLH, metres and m/s.
% The browser handles rendering/playback; no graphics toolbox or server needed.

if nargin < 3
    target = [50; 50; 0];
end
validateattributes(time, {'numeric'}, {'real','finite','vector','nonempty'});
validateattributes(state, {'numeric'}, {'real','finite','2d','ncols',6});
validateattributes(target, {'numeric'}, {'real','finite','vector','numel',3});
time = double(time(:));
assert(numel(time) >= 2 && size(state,1) == numel(time), ...
    'Provide at least two time samples and one state row per time sample.');
assert(all(diff(time) > 0), 'Time must be strictly increasing.');

% Limit browser data size for long simulations. This affects display only.
indices = unique(round(linspace(1, numel(time), min(numel(time), 12000))));
data.time = time(indices);
data.state = double(state(indices,:));
data.target = double(target(:).');
data.demo = false;

% Homework leader orbit. Keep a and mu consistent with CW_Matrix.m.
data.orbit = struct('a',6898e3, 'mu',3.986e14, ...
    'inclinationDeg',23, 'raanDeg',10, 'argPeriapsisDeg',90, ...
    'trueAnomaly0Deg',0, 'earthRadius',6378.1e3);

% Reconstruct ECI positions in the renderer:
% Q = Rz(RAAN)*Rx(inclination)*Rz(argument of periapsis).
% theta = nu0 + sqrt(mu/a^3)*(time-time(1)).
% eR = Q*[cos(theta);sin(theta);0], eT = Q*[-sin(theta);cos(theta);0].
% eN = Q*[0;0;1]; rLeader = a*eR; rFollower = rLeader + [eR eT eN]*rho.
% Bodies/panels are illustrative; this project does not simulate attitude.

folder = fileparts(mfilename('fullpath'));
template = fileread(fullfile(folder, 'Orbit_Renderer_Template.html'));
token = '__PROJECT4_DATA_JSON__';
assert(numel(strfind(template, token)) == 1, ...
    'The renderer template must contain exactly one data placeholder.');
html = strrep(template, token, jsonencode(data));
htmlFile = fullfile(folder, 'Orbit_Renderer.html');
[fid, message] = fopen(htmlFile, 'w', 'n', 'UTF-8');
assert(fid ~= -1, 'Cannot write the renderer: %s', message);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', html);
clear cleanup; % Close the file before the browser opens it.
fprintf('3-D renderer saved: %s\n', htmlFile);
status = web(htmlFile, '-browser');
if status ~= 0
    fprintf('Open Orbit_Renderer.html from this folder in Chrome or Edge.\n');
end
end
