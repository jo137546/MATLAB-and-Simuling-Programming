

%%% initialize %%%

% Clear workspace to avoid variable overwrites
clear all;

% Initialize
h = 15000; % altitude in meters
Ma = 1.8; % Mach number
Ta = 1500; % turbine inlet max temperature (K)
LHV = 43124000; % LHV (J/kg)
fst = 0.06; % stoichiometric fuel-air ratio

nud = 0.9; % diffuser efficiency
nub = 0.98; % combustion chamber efficiency
nut = 0.92; % turbine efficiency
nuc = 0.9; % compressor efficiency
nun = 0.98; % nozzle efficiency

rb = 0.97; % pressure ratio across combustion chamber

gam1 = 1.4; % heat capacity ratio (to burner)
gam2 = 1.3; % heat capacity ratio (rest of engine)
R = 287; % gas constant (J/kg·K)

% Ambient conditions at 15,000 m
Tamb = 273.15 - 56.7; % K
Pamb = 12100; % Pa

% Afterburner constants
T06_ab = 2000; eta_ab = 0.95; rab = 0.9;

% RC sweep setup
rc_vals = 2:1:60;
n = length(rc_vals);

% Preallocate w/ no afterburner
Fxs_vals = zeros(1,n);
TSFC_vals = zeros(1,n);
nuth_vals = zeros(1,n);
nuo_vals = zeros(1,n);
nup_vals = zeros(1,n);
A7A6_vals = zeros(1, n); % Nozzle area ratio A7/A6

%preallocate w/ afterburner
Fxs_ab_vals = zeros(1,n);
TSFC_ab_vals = zeros(1,n);
nuth_ab_vals = zeros(1,n);
nuo_ab_vals = zeros(1,n);
nup_ab_vals = zeros(1,n);
A7A6_ab_vals = zeros(1, n); % Nozzle area ratio A7/A6

%%% begin calculations %%%
for i = 1:n
    rc = rc_vals(i);

    % Freestream (Station 0)
    T0amb = Tamb * (1 + ((gam1 - 1) / 2) * Ma^2); % Stagnation temperature at inlet
    P0amb = Pamb * ((1 + ((gam1 - 1) / 2) * Ma^2) ^ (gam1 / (gam1 - 1))); % Stagnation pressure at inlet
    u = Ma * sqrt(gam1 * R * Tamb); % Engine velocity through air
    
    % Diffuser (Station 0 to 2)
    T02 = T0amb; % Adiabatic, so T02 = T0amb
    T02s = (nud * (T0amb - Tamb)) + Tamb;
    P02 = Pamb * (T02s / Tamb) ^ (gam1 / (gam1 - 1));
    
    % Compressor (Station 2 to 3)
    P03 = rc * P02;
    gam1 = 1.4; % Redefine to be sure
    exponent1 = (gam1 - 1) / gam1;
    term = (P03 / P02) ^ exponent1;
    T03s = T02 * term;
    T03 = ((T03s - T02) / nuc) + T02;
    
    Cp1 = R / (1 - (1 / gam1));
    Wc = Cp1 * (T03 - T02); % finding work generated in compressor
    
    % Combustor/burner (Station 3 to 4)
    T04 = Ta;
    Cp2 = R / (1 - (1 / gam2));
    fb = (Cp2 * (T04 - T03)) / (nub * LHV - Cp2 * T04);
    
    fb_min = 0.003;

    if fb > fst || fb < fb_min

        warning('Skipping unphysical case: fb = %.5f', fb);
        Fxs_vals(i) = NaN;
        TSFC_vals(i) = NaN;
        nuth_vals(i) = NaN;
        nuo_vals(i) = NaN;
        nup_vals(i) = NaN;
        A7A6_vals(i) = NaN;
        continue;

    end
    
    P04 = rb * P03;
    
    % Turbine (station 4 to 5)
    exponent2 = (gam2 / (gam2 - 1));
    Wtout = Wc;
    T04_minus_T05 = Wtout / (Cp2 * (1 + fb));
    T05 = T04 - T04_minus_T05;
    T05s = T04 - ((T04 - T05) / nut);
    P05 = P04 * ((T05s / T04) ^ exponent2);
    
        % nozzle NO AB (station 5 to 7)
        T06 = T05;
        T07 = T06;
        P06 = P05; %negligible losses
        P7 = Pamb;
        
        exponent3 = ((gam2 - 1) / gam2);
        
        T7s = T06 / ((P06 / P7) ^ exponent3);
        T7 = T06 - (nun * (T06 - T7s));
        
        M7 = sqrt(((T07 / T7) - 1) / ((gam2 - 1) / 2 )); %exit mach
        u7 = M7 * sqrt(gam1 * R * T7); %velocity at nozzle exit
        
        %finding TSFC
        ue = sqrt(2 * nun * Cp2 * (T06 - T7s)); % Exit velocity (ideal expansion)
        Fxs = ue - u; % Specific thrust (N·s/kg)
        TSFC = fb / Fxs;
        
        %finding thermal efficiency
        nuth = (((1 + fb) * ((ue^2)/2)) - (u^2 /2)) / (fb * LHV);
        
        %finding overall efficiency
        nuo = (Fxs * u) / (fb * LHV);
        
        %finding propulsive efficiency
        nup = nuo / nuth;
    
        gamma_n = gam2; % gas properties after turbine
        area_ratio = (1 / M7) * ((2 / (gamma_n + 1)) * (1 + ((gamma_n - 1)/2) * M7^2)) ^ ((gamma_n + 1) / (2 * (gamma_n - 1)));
      
        % Store results for no AB
        Fxs_vals(i) = Fxs;
        TSFC_vals(i) = TSFC;
        nuth_vals(i) = nuth;
        nuo_vals(i) = nuo;
        nup_vals(i) = nup;
        A7A6_vals(i) = area_ratio;
    
        % nozzle with AB
        f_ab = (Cp2 * (T06_ab - T05)) / (eta_ab * LHV - Cp2 * T06_ab);
        f_total = fb + f_ab;
    
        if f_ab < 0 || f_total > fst
            continue;
        end

        T06 = T06_ab;
        P06 = rab * P05;
    
        % Nozzle (afterburner path)
        T7s_ab = T06 / (P06 / Pamb)^((gam2 - 1)/gam2);
        T7_ab = T06 - nun * (T06 - T7s_ab);
        ue_ab = sqrt(2 * nun * Cp2 * (T06 - T7s_ab));
        M7_ab = sqrt(((T06 / T7_ab) - 1) * 2 / (gam2 - 1));
        area_ratio_ab = (1 / M7_ab) * ((2 / (gam2 + 1)) * (1 + ((gam2 - 1)/2) * M7_ab^2))^((gam2 + 1) / (2 * (gam2 - 1)));
        Fxs_ab = ue_ab - u;
        TSFC_ab = f_total / Fxs_ab;
        nuth_ab = (((1 + f_total) * ue_ab^2 / 2) - u^2 / 2) / (f_total * LHV);
        nuo_ab = (Fxs_ab * u) / (f_total * LHV); nup_ab = nuo_ab / nuth_ab;
    
        % Store AB results
        Fxs_ab_vals(i) = Fxs_ab; TSFC_ab_vals(i) = TSFC_ab;
        nuth_ab_vals(i) = nuth_ab; nuo_ab_vals(i) = nuo_ab;
        nup_ab_vals(i) = nup_ab; A7A6_ab_vals(i) = area_ratio_ab;
end

% % for debugging
% fprintf('T02 = %.2f K\n', T02);
% fprintf('P02 = %.2f Pa\n', P02);
% fprintf('T03s = %.2f K\n', T03s);
% fprintf('T03 = %.2f K\n', T03);
% fprintf('Cp1 = %.2f J/kg·K\n', Cp1);
% fprintf('Wc = %.4f J/kg\n', Wc);
% fprintf('Cp2 = %.2f J/kg·K\n', Cp2);
% fprintf('fb = %.4f\n', fb);
% fprintf('P04 = %.4f Pa\n', P04);
% fprintf('T05 = %.4f K\n', T05);
% fprintf('T05s = %.4f K\n', T05s);
% fprintf('P05 = %.4f Pa\n', P05);
% fprintf('T7s = %.4f K\n', T7s);
% fprintf('T7 = %.4f K\n', T7);
% fprintf('M7 = %.4f\n', M7);
% fprintf('u7 = %.4f m/s\n', u7);
% fprintf('ue = %.4f m/s\n', ue);
% fprintf('Specific thrust = %.4f N·s/kg\n', Fxs);
% fprintf('TSFC = %.10f  kg/N·s \n', TSFC);
% fprintf('nuth = %.4f\n', nuth);
% fprintf('nuo = %.4f\n', nuo);
% fprintf('nup = %.4f\n', nup); 

%%% plotting %%%

% 1. Specific Thrust
figure;
plot(rc_vals, Fxs_vals, 'b', rc_vals, Fxs_ab_vals, 'r--', 'LineWidth', 1.5);
xlabel('Compressor Pressure Ratio (r_c)'); ylabel('Specific Thrust (N·s/kg)');
legend('No Afterburner','With Afterburner'); title('Specific Thrust vs r_c');
grid on;

% 2. TSFC
figure;
plot(rc_vals, TSFC_vals*1e6, 'b', rc_vals, TSFC_ab_vals*1e6, 'r--', 'LineWidth', 1.5);
xlabel('Compressor Pressure Ratio (r_c)'); ylabel('TSFC (mg/N·s)');
legend('No Afterburner','With Afterburner'); title('TSFC vs r_c');
grid on;

% 3. Efficiencies
figure;
plot(rc_vals, nuth_vals, 'r', rc_vals, nup_vals, 'g', rc_vals, nuo_vals, 'b', ...
     rc_vals, nuth_ab_vals, 'r--', rc_vals, nup_ab_vals, 'g--', rc_vals, nuo_ab_vals, 'b--', 'LineWidth', 1.5);
xlabel('r_c'); ylabel('Efficiency'); title('Efficiencies vs r_c');
legend('\eta_{th}','\eta_p','\eta_o','\eta_{th,AB}','\eta_{p,AB}','\eta_{o,AB}');
ylim([0 1]); grid on;

% 4. Nozzle Area Ratio
figure;
plot(rc_vals, A7A6_vals, 'b', rc_vals, A7A6_ab_vals, 'r--', 'LineWidth', 1.5);
xlabel('Compressor Pressure Ratio (r_c)'); ylabel('Nozzle Area Ratio (A_7 / A_6)');
legend('No Afterburner','With Afterburner'); title('Nozzle Area Ratio vs r_c');
grid on;