function mdl = build_ss_simulink(A,B,C,D, varargin)
% BUILD_SS_SIMULINK  Create a Simulink model from A,B,C,D.
% Usage:
%   mdl = build_ss_simulink(A,B,C,D)               % continuous
%   mdl = build_ss_simulink(A,B,C,D,'Ts',0.01)     % discrete with Ts
%   mdl = build_ss_simulink(A,B,C,D,'Name','myModel','Harness',true)
%
% Options:
%   'Ts'      : sample time. If provided and >0, uses Discrete State-Space.
%   'Name'    : model name (default 'ss_autogen')
%   'Harness' : true adds Step→System→Scope for quick testing (default true)

opts.Ts = [];
opts.Name = 'ss_autogen';
opts.Harness = true;
for k = 1:2:numel(varargin)
    opts.(varargin{k}) = varargin{k+1};
end

% Basic checks
[nx,nu] = size(B); 
[ny,~]  = size(C);
assert(all(size(A) == [nx nx]), 'A must be nx-by-nx.');
assert(all(size(D) == [ny nu]), 'D must be ny-by-nu.');

% Make / reset model
mdl = opts.Name;
if bdIsLoaded(mdl); close_system(mdl,0); end
new_system(mdl); open_system(mdl);

% Choose block type
if ~isempty(opts.Ts) && opts.Ts > 0
    blkType = 'simulink/Discrete/Discrete State-Space';
else
    blkType = 'simulink/Continuous/State-Space';
end

% Add the State-Space block to the model
ssBlk = add_block(blkType, [mdl '/Plant'], 'Position', [200 80 300 200]);

for i = 1:nu
    add_block('simulink/Sources/In1', [mdl sprintf('/u%d', i)], 'Position', [40 80 + (i-1)*40 70 100 + (i-1)*40]);
end

for i = 1:ny
    add_block('simulink/Sinks/Out1', [mdl sprintf('/y%d', i)], 'Position', [420 80 + (i-1)*40 450 100 + (i-1)*40]);
end

% Set A,B,C,D (and Ts for discrete)
set_param(ssBlk, ...
    'A', mat2str(A), ...
    'B', mat2str(B), ...
    'C', mat2str(C), ...
    'D', mat2str(D));

if ~isempty(opts.Ts) && opts.Ts > 0
    set_param(ssBlk,'SampleTime', num2str(opts.Ts));
end

% Wire main path
add_line(mdl,'u/1','Plant/1','autorouting','on');
add_line(mdl,'Plant/1','y/1','autorouting','on');

% Optional harness: Step + Scope for quick testing
if opts.Harness
    stepBlk  = add_block('simulink/Sources/Step', [mdl '/Step'], 'Position',[40 20 70 40]);
    scopeBlk = add_block('simulink/Sinks/Scope', [mdl '/Scope'], 'Position',[420 20 450 40]);

    % Set block parameters to avoid warnings
    set_param(stepBlk, 'Time', '0', 'Before', '1', 'After', '1');
    set_param(scopeBlk, 'NumInputPorts', '1');

    % For MIMO, Step is scalar; it drives port 1. You can add more sources if needed.
    add_line(mdl,'Step/1','u/1','autorouting','on');
    add_line(mdl,'Plant/1','Scope/1','autorouting','on');
end

% Tidy up layout
try 
    Simulink.BlockDiagram.arrangeSystem(mdl); 
catch ME
    warning('Failed to arrange the system: %s', E.message);
    % Optionally, you can add a command to display the model for debugging
    open_system(mdl);  % Open the model to visually inspect the layout
end

% Save (optional)
save_system(mdl);
open_system(mdl);  % bring to front
end