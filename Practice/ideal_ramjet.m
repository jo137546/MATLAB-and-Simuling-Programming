% Ask user for unit system
unit_system = input('Enter unit system (1 for Imperial, 2 for Metric): ');
if unit_system ~= 1 && unit_system ~= 2
    error('Invalid unit system selection. Please run again and select 1 or 2.');
end

% Define parameters based on unit system
if unit_system == 1    % Imperial
    g = 32.2;         % ft/s^2
    R = 1718;         % ft·lbf/(slug·°R)
    unit_str = {'PSI', 'ft'};
    % Conversion factors
    m_to_ft = 3.28084; % m to ft
    kJkg_to_Btulb = 430; % 1 kJ/kg = 430 Btu/lb (approximate for jet fuel)
    K_to_R = 1.8;      % Kelvin to Rankine
else    % Metric
    g = 9.81;         % m/s^2
    R = 287;          % J/(kg·K)
    unit_str = {'Pa', 'm'};
    % Conversion factors (for reference, though mostly not needed in Metric)
    m_to_ft = 1;       % No conversion needed
    kJkg_to_Btulb = 1; % No conversion needed
    K_to_R = 1;        % No conversion needed
end

% Altitude (given as 10,000 m, convert if Imperial)
altitude = 10000; % m (as given in the problem)
if unit_system == 1
    altitude = altitude * m_to_ft; % Convert to ft
end
fprintf('Altitude: %.2f %s\n', altitude, unit_str{2});

% Specific heat of combustion (given as 43,000 kJ/kg, convert if Imperial)
if unit_system == 1
    delh = 43000 * kJkg_to_Btulb; % Convert to Btu/lb
    fprintf('Heat of combustion: %.2f Btu/lb (converted from 43,000 kJ/kg)\n', delh);
else
    delh = 43000e3; % J/kg (43,000 kJ/kg as given)
    fprintf('Heat of combustion: %.2f J/kg\n', delh);
end

% Fuel-to-air ratio (given as stoichiometric 0.06)
f = 0.06; % Stoichiometric fuel-to-air ratio
fprintf('Using stoichiometric fuel-to-air ratio: %.3f\n', f);

% Ramjet max temperature (given as 2600 K, convert if Imperial)
T04 = 2600; % K (maximum temperature in combustor)
if unit_system == 1
    T04_imperial = T04 * K_to_R; % Convert to °R
else
    T04_imperial = T04; % Keep in K for Metric
end
% Fix the ternary operator syntax (Line 53)
if unit_system == 1
    temp_unit = '°R';
else
    temp_unit = 'K';
end
fprintf('Using maximum ramjet temperature: %.1f %s\n', T04_imperial, temp_unit);

% Specific heat ratio (given as 1.4)
gamma = 1.4;

% Specific heat at constant pressure (calculated using given equation)
Cp = (gamma / (gamma - 1)) * R;

% Atmospheric conditions at 10,000 m (from ISA model, Appendix III)
if unit_system == 2
    T0 = 223.15; % K at 10,000 m (ISA standard)
    P0 = 26436;  % Pa at 10,000 m (ISA standard)
    rho0 = 0.660; % kg/m^3 at 10,000 m (ISA standard)
else
    T0 = 223.15 * K_to_R; % Convert to °R (Rankine)
    P0 = 26436 * 0.000145038; % Convert Pa to PSI
    rho0 = 0.660 * 0.06242796; % Convert kg/m^3 to slug/ft^3
end
a0 = sqrt(gamma * R * T0); % Speed of sound at altitude
% Fix the ternary operator syntax for speed of sound
if unit_system == 1
    speed_unit = 'ft/s';
else
    speed_unit = 'm/s';
end
fprintf('Speed of sound at altitude: %.2f %s\n', a0, speed_unit);

% Flight Mach number range (1 to 6, with step size for tabulation)
M_flight = 1:1:6; % Step size of 1 for tabulation (1, 2, 3, 4, 5, 6)
M_flight_plot = 1:0.1:6; % Finer step for plotting
n = length(M_flight);
n_plot = length(M_flight_plot);

% Initialize arrays for results (for tabulation)
specific_thrust = zeros(1, n);
TSFC = zeros(1, n);
T04_calc = zeros(1, n);
A_exit_A_throat = zeros(1, n);
eta_th = zeros(1, n);
eta_p = zeros(1, n);
eta_0 = zeros(1, n);

% Initialize arrays for plotting
specific_thrust_plot = zeros(1, n_plot);
TSFC_plot = zeros(1, n_plot);
T04_calc_plot = zeros(1, n_plot);
A_exit_A_throat_plot = zeros(1, n_plot);
eta_th_plot = zeros(1, n_plot);
eta_p_plot = zeros(1, n_plot);
eta_0_plot = zeros(1, n_plot);

% Loop over Mach numbers for tabulation
for i = 1:n
    M0 = M_flight(i);
    u0 = M0 * a0; % Freestream velocity

    % Stagnation conditions at inlet (station 0 to 2, assuming ideal inlet)
    T02 = T0 * (1 + (gamma - 1)/2 * M0^2);
    P02 = P0 * (1 + (gamma - 1)/2 * M0^2)^(gamma/(gamma-1));

    % Energy balance in combustor (station 2 to 4)
    % Cp * (T04 - T02) = f * delh
    T04_calc(i) = T02 + (f * delh) / Cp;

    % Check if T04 exceeds maximum temperature
    if unit_system == 1
        T04_limit = T04 * K_to_R; % Convert max temp to °R for comparison
    else
        T04_limit = T04; % Keep in K
    end
    if T04_calc(i) > T04_limit
        % Adjust fuel-to-air ratio to meet T04 limit
        f_adjusted = (T04_limit - T02) * Cp / delh;
    else
        f_adjusted = f;
    end

    % Exit temperature (station 4 to 9, assuming ideal nozzle)
    % For simplicity, assume T9 = T0 (fully expanded nozzle)
    % Use energy equation to find exit velocity
    u9 = sqrt(2 * Cp * (T04_limit - T0));

    % Specific thrust: (u9 - u0) / g
    specific_thrust(i) = (u9 - u0) / g;

    % TSFC: f / specific_thrust
    TSFC(i) = f_adjusted / specific_thrust(i);
    if unit_system == 1
        TSFC(i) = TSFC(i) * (1 / 3600); % Convert to lb/lbf·s (approximate)
    else
        TSFC(i) = TSFC(i) * (1 / 1000); % Convert to kg/N·s (approximate)
    end

    % Area ratio (A_exit/A_throat) - simplified
    M9 = u9 / sqrt(gamma * R * T0);
    term1 = 1 / M9;
    term2 = 2 / (gamma + 1);
    term3 = 1 + (gamma - 1) / 2 * M9^2;
    exponent = (gamma + 1) / (2 * (gamma - 1));
    A_exit_A_throat(i) = term1 * (term2 * term3)^exponent;

    % Efficiencies
    eta_th(i) = 1 - (T0 / T04_limit); % Thermal efficiency
    eta_p(i) = 2 * u0 / (u9 + u0); % Propulsive efficiency
    eta_0(i) = eta_th(i) * eta_p(i); % Overall efficiency
end

% Convert T04_calc to appropriate units for tabulation
if unit_system == 1
    T04_calc = T04_calc * K_to_R; % Convert to °R
end

% Display tabulated results
fprintf('\nTabulated Results for Mach Number 1 to 6:\n');
if unit_system == 2
    fprintf('Mach\tSpec. Thrust (N·s/kg)\tTSFC (kg/N·s)\tT04 (K)\tA_exit/A_throat\teta_th\teta_p\teta_0\n');
else
    fprintf('Mach\tSpec. Thrust (lbf·s/slug)\tTSFC (lb/lbf·s)\tT04 (°R)\tA_exit/A_throat\teta_th\teta_p\teta_0\n');
end
for i = 1:n
    fprintf('%.1f\t%.2f\t\t%.6f\t\t%.2f\t%.2f\t\t%.3f\t%.3f\t%.3f\n', ...
        M_flight(i), specific_thrust(i), TSFC(i), T04_calc(i), A_exit_A_throat(i), ...
        eta_th(i), eta_p(i), eta_0(i));
end

% Loop over Mach numbers for plotting (finer resolution)
for i = 1:n_plot
    M0 = M_flight_plot(i);
    u0 = M0 * a0; % Freestream velocity

    % Stagnation conditions at inlet
    T02 = T0 * (1 + (gamma - 1)/2 * M0^2);
    P02 = P0 * (1 + (gamma - 1)/2 * M0^2)^(gamma/(gamma-1));

    % Energy balance in combustor
    T04_calc_plot(i) = T02 + (f * delh) / Cp;

    % Check if T04 exceeds maximum temperature
    if T04_calc_plot(i) > T04_limit
        f_adjusted = (T04_limit - T02) * Cp / delh;
    else
        f_adjusted = f;
    end

    % Exit velocity
    u9 = sqrt(2 * Cp * (T04_limit - T0));

    % Specific thrust
    specific_thrust_plot(i) = (u9 - u0) / g;

    % TSFC
    TSFC_plot(i) = f_adjusted / specific_thrust_plot(i);
    if unit_system == 1
        TSFC_plot(i) = TSFC_plot(i) * (1 / 3600);
    else
        TSFC_plot(i) = TSFC_plot(i) * (1 / 1000);
    end

    % Area ratio
    M9 = u9 / sqrt(gamma * R * T0);
    term1 = 1 / M9;
    term2 = 2 / (gamma + 1);
    term3 = 1 + (gamma - 1) / 2 * M9^2;
    exponent = (gamma + 1) / (2 * (gamma - 1));
    A_exit_A_throat_plot(i) = term1 * (term2 * term3)^exponent;

    % Efficiencies
    eta_th_plot(i) = 1 - (T0 / T04_limit);
    eta_p_plot(i) = 2 * u0 / (u9 + u0);
    eta_0_plot(i) = eta_th_plot(i) * eta_p_plot(i);
end

% Convert T04_calc_plot to appropriate units for plotting
if unit_system == 1
    T04_calc_plot = T04_calc_plot * K_to_R; % Convert to °R
end

% Plotting
% a) Specific thrust vs. M_flight
figure;
plot(M_flight_plot, specific_thrust_plot, 'b-', 'LineWidth', 2);
xlabel('Flight Mach Number (M_{flight})');
if unit_system == 2
    ylabel('Specific Thrust (N·s/kg)');
else
    ylabel('Specific Thrust (lbf·s/slug)');
end
title('Specific Thrust vs. Flight Mach Number');
grid on;

% b) TSFC vs. M_flight
figure;
plot(M_flight_plot, TSFC_plot, 'r-', 'LineWidth', 2);
xlabel('Flight Mach Number (M_{flight})');
if unit_system == 2
    ylabel('TSFC (kg/N·s)');
else
    ylabel('TSFC (lb/lbf·s)');
end
title('TSFC vs. Flight Mach Number');
grid on;

% c) T04 vs. M_flight
figure;
plot(M_flight_plot, T04_calc_plot, 'g-', 'LineWidth', 2);
xlabel('Flight Mach Number (M_{flight})');
if unit_system == 2
    ylabel('Combustor Exit Temperature T_{04} (K)');
else
    ylabel('Combustor Exit Temperature T_{04} (°R)');
end
title('T_{04} vs. Flight Mach Number');
grid on;

% d) A_exit/A_throat vs. M_flight
figure;
plot(M_flight_plot, A_exit_A_throat_plot, 'm-', 'LineWidth', 2);
xlabel('Flight Mach Number (M_{flight})');
ylabel('A_{exit}/A_{throat}');
title('Exit-to-Throat Area Ratio vs. Flight Mach Number');
grid on;

% e) Efficiencies vs. M_flight
figure;
plot(M_flight_plot, eta_th_plot, 'b-', 'LineWidth', 2, 'DisplayName', '\eta_{th}');
hold on;
plot(M_flight_plot, eta_p_plot, 'r-', 'LineWidth', 2, 'DisplayName', '\eta_{p}');
plot(M_flight_plot, eta_0_plot, 'g-', 'LineWidth', 2, 'DisplayName', '\eta_{0}');
xlabel('Flight Mach Number (M_{flight})');
ylabel('Efficiency');
title('Efficiencies vs. Flight Mach Number');
legend('show');
grid on;