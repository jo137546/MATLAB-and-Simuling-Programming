% ===== file: build_ss_autowire.m =====
% Requires A,B,C,D,Ts in base workspace (run init_ss first).
function mdl = build_ss_autowire(modelName)
if nargin<1, modelName = 'ss_autogen'; end

% Create/overwrite model
if bdIsLoaded(modelName), close_system(modelName, 0); end
new_system(modelName); open_system(modelName);

try
    % Add blocks
    stepBlk = add_block('simulink/Sources/Step', [modelName '/Step'], 'Position',[40 80 70 110]);
    if isempty(evalin('base','Ts'))
        ssLib = 'simulink/Continuous/State-Space';
    else
        ssLib = 'simulink/Discrete/Discrete State-Space';
    end
    ssBlk   = add_block(ssLib, [modelName '/Plant'], 'Position',[180 65 340 125]);
    scopeBlk= add_block('simulink/Sinks/Scope', [modelName '/Scope'], 'Position',[410 80 440 110]);

    % Point block params at base-workspace variables
    set_param(ssBlk, 'A','A', 'B','B', 'C','C', 'D','D');
    if isempty(evalin('base','Ts'))
        set_param(ssBlk, 'SampleTime','');   % continuous
    else
        set_param(ssBlk, 'SampleTime','Ts'); % discrete uses Ts
    end

    % Wire: Step -> Plant -> Scope
    add_line(modelName, 'Step/1', 'Plant/1','autorouting','on');
    add_line(modelName, 'Plant/1', 'Scope/1','autorouting','on');

    % Nice default step (optional)
    set_param(stepBlk, 'Time','0', 'Before','0', 'After','1');

    save_system(modelName);
    mdl = modelName;
    fprintf('Model "%s" ready. Press Run.\n', modelName);

catch ME
    close_system(modelName, 0);
    rethrow(ME);
end
end