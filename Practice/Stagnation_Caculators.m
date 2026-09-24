%%%%%% Initialization %%%%%%

% Ask user for unit system
unit_system = input('Enter unit system (1 for Imperial, 2 for Metric): ');
if unit_system ~= 1 && unit_system ~= 2
    error('Invalid unit system selection. Please run again and select 1 or 2.');
end

% Define parameters based on unit system
if unit_system == 1    % Imperial
    params = {'Velocity (ft/s)', 'u', 'Height (ft)', 'h'};
    g = 32.2;         % ft/s^2
    h_max_tropo = 36152;  % ft
    h_max_lower_strato = 82345; % ft
    P0 = 2116;        % psf
    T0 = 59;          % °F
    lapse_rate = 0.00356;  % °F/ft
    R = 1718;         % ft·lbf/(slug·°R)
    unit_str = {'PSI', 'ft'};
else    % Metric
    params = {'Velocity (m/s)', 'u', 'Height (m)', 'h'};
    g = 9.81;         % m/s^2
    h_max_tropo = 11000;  % m
    h_max_lower_strato = 25100; % m
    P0 = 101325;      % Pa
    T0 = 15;          % °C
    lapse_rate = 0.0065;   % °C/m
    R = 287;          % J/(kg·K)
    unit_str = {'Pa', 'm'};
end

% Get input with error checking
values = zeros(1,2);
for i = 1:2:length(params)
    while true
        prompt = sprintf('Enter %s: ', params{i});
        val = input(prompt);
        if val < 0
            fprintf('Error: Value must be non-negative. Try again.\n');
        else
            values((i+1)/2) = val;
            break;
        end
    end
end

u = values(1);    % velocity
h = values(2);    % height

%%%%%% Pressure Calculations %%%%%%

% Convert to Imperial for calculations (since equations are in Imperial)
h_imp = h * (unit_system == 2) * 3.28084 + h * (unit_system == 1);  % m to ft if metric
u_imp = u * (unit_system == 2) * 3.28084 + u * (unit_system == 1);  % m/s to ft/s if metric

% Calculate temperature and pressure based on height (using Imperial equations)
if h_imp < 36152    % Troposphere
    T_imp = 59 - 0.00356 * h_imp;
    T_abs = T_imp + 459.7;
    P_static = 2116 * ((T_imp + 459.7) / 518.7)^5.256;  % Corrected reference temperature
elseif h_imp < 82345    % Lower Stratosphere
    T_imp = -70;
    T_abs = T_imp + 459.7;
    P_static = 473.1 * exp(1.73 - 0.000048 * h_imp);
else    % Upper Stratosphere
    T_imp = -205.05 + 0.00164 * h_imp;
    T_abs = T_imp + 459.7;
    P_static = 51.97 * ((T_imp + 459.7) / 389.98)^-11.388;
end

% Calculate density
rho_imp = P_static / (1718 * (T_imp + 459.7));

% Convert back to selected unit system for density
rho = rho_imp * (unit_system == 1) + rho_imp * (unit_system == 2) * 515.379;  % slug/ft³ to kg/m³ if metric
P_static = P_static * (unit_system == 1) + P_static * (unit_system == 2) * 47.88026;  % psf to Pa if metric
dynamic_pressure = 0.5 * rho * u^2;
P_total = P_static + dynamic_pressure;

% Convert units for display in selected system
if unit_system == 1
    P_static_display = P_static / 144;    % psf to psi
    P_dynamic_display = dynamic_pressure / 144;
    P_total_display = P_total / 144;
    unit_display = 'PSI';
else
    P_static_display = P_static;    % Already in Pa
    P_dynamic_display = dynamic_pressure;
    P_total_display = P_total;
    unit_display = 'Pa';
end

%%%%%% Table Output %%%%%%

% Define test points
if unit_system == 1
    heights = [0, 10000, 20000, 30000, 40000];
    velocities = [0, 100, 200, 300, 400];
    temp_unit = '°F';
    dens_unit = 'slug/ft³';
else
    heights = [0, 3000, 6000, 9000, 12000];
    velocities = [0, 30, 60, 90, 120];
    temp_unit = '°C';
    dens_unit = 'kg/m³';
end

fprintf('\nAtmospheric Properties Table:\n');
fprintf('----------------------------------------------------------------------------------------\n');
fprintf('Height (%s)  Vel (%s)  Temp (%s)  Static P (%s)  Density (%s)  Dynamic P (%s)\n', ...
    unit_str{2}, unit_str{2}(1), temp_unit, unit_display, dens_unit, unit_display);
fprintf('----------------------------------------------------------------------------------------\n');

for h_test = heights
    for u_test = velocities
        % Convert to Imperial for calculations
        h_test_imp = h_test * (unit_system == 2) * 3.28084 + h_test * (unit_system == 1);
        u_test_imp = u_test * (unit_system == 2) * 3.28084 + u_test * (unit_system == 1);
        
        % Recalculate for each point
        if h_test_imp < 36152
            T = 59 - 0.00356 * h_test_imp;
            P_s = 2116 * ((T + 459.7) / 518.7)^5.256;  % Corrected reference temperature
        elseif h_test_imp < 82345
            T = -70;
            P_s = 473.1 * exp(1.73 - 0.000048 * h_test_imp);
        else
            T = -205.05 + 0.00164 * h_test_imp;
            P_s = 51.97 * ((T + 459.7) / 389.98)^-11.388;
        end
        
        rho_test = P_s / (1718 * (T + 459.7));
        rho_test = rho_test * (unit_system == 1) + rho_test * (unit_system == 2) * 515.379;
        P_s = P_s * (unit_system == 1) + P_s * (unit_system == 2) * 47.88026;
        P_d = 0.5 * rho_test * u_test^2;
        
        % Convert temperature for display
        T = T * (unit_system == 1) + (T - 32) * 5/9 * (unit_system == 2);
        
        % Convert units for display
        if unit_system == 1
            P_s_display = P_s / 144;
            P_d_display = P_d / 144;
        else
            P_s_display = P_s;
            P_d_display = P_d;
        end
        
        fprintf('%8d %8d %10.1f %12.2f %14.5f %12.2f\n', ...
            h_test, u_test, T, P_s_display, rho_test, P_d_display);
    end
end

% Single point output with both Imperial and Metric
fprintf('\nSingle Point Calculation:\n');
fprintf('Input: h = %d %s, u = %d %s\n', h, unit_str{2}, u, unit_str{2}(1));

% Imperial units (already calculated)
fprintf('\nImperial Units:\n');
fprintf('Height: %.1f ft\n', h_imp);
fprintf('Velocity: %.1f ft/s\n', u_imp);
fprintf('Static Temperature: %.1f °F\n', T_imp);
fprintf('Static Pressure: %.2f PSI\n', P_static / 144);
fprintf('Density: %.5f slug/ft³\n', rho_imp);
fprintf('Dynamic Pressure: %.2f PSI\n', dynamic_pressure / 144);
fprintf('Total Pressure: %.2f PSI\n', P_total / 144);

% Metric units
h_met = h_imp * 0.3048;    % ft to m
u_met = u_imp * 0.3048;    % ft/s to m/s
T_met = (T_imp - 32) * 5/9;  % °F to °C
P_static_met = (P_static / 144) * 6894.76;  % psi to Pa
P_dynamic_met = (dynamic_pressure / 144) * 6894.76;
P_total_met = (P_total / 144) * 6894.76;
rho_met = rho_imp * 515.379;  % slug/ft³ to kg/m³

fprintf('\nMetric Units:\n');
fprintf('Height: %.1f m\n', h_met);
fprintf('Velocity: %.1f m/s\n', u_met);
fprintf('Static Temperature: %.1f °C\n', T_met);
fprintf('Static Pressure: %.2f Pa\n', P_static_met);
fprintf('Density: %.5f kg/m³\n', rho_met);
fprintf('Dynamic Pressure: %.2f Pa\n', P_dynamic_met);
fprintf('Total Pressure: %.2f Pa\n', P_total_met);