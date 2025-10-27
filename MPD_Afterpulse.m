%% MicroPulse DIAL Afterpulse Error Simulation 

clear; close all; clc;


%% Initialization and Unit Conversion

% System Parameters
tau_pulse = 200e-9;         % Pulse width (s)
R_start = 100;              % Start range for retrieval (m)
R_end = 4000;               % End range for retrieval (m)
dR = 7.5;                   % Range resolution (m/bin)
P_laser = 50e-3;            % Avg laser power (W) 
sigma_on_cm2 = 0.99e-23;    % Online cross-section (cm^2/molecule): [use for 20 g/m^3 wv surf] 
%sigma_on_cm2 = 2.5e-23;    % Online cross-section (cm^2/molecule): [use for 7.5 g/m^3 wv surf] 
sigma_off_cm2 = 7.0e-25;    % Offline cross-section (cm^2/molecule)

% Detector/Afterpulse Parameters
N_dark = 0;                  % Constant Dark Count Rate (counts/bin)
lambda_AP = 50;             % Afterpulse decay length (m)
AP_fraction = 1e-13;         % Afterpulse peak as fraction of C_sys

% Geometric Overlap Parameters
R_full_overlap = 1500;       % Range at which O(R) approaches 1 (m)
O_scale = 100;               % Scaling factor for the exponential rise 

% Simple Atmosphere Parameters 
N_solar = 0;              % Solar background Count Rate (counts/bin)
alpha_non_abs = 1e-8;     % Total aerosol and molecular (non-absorbing) extinction (m^-1)
wv_mass_surf = 20;        % Water Vapor Density Profile
%wv_mass_surf = 7.5;      % Water Vapor Density Profile
wv_mass_true = @(R) wv_mass_surf .* exp(-R/2000); % Mass density profile (g/m^3)
beta_a = 1e-7;            % general backscatter coefficient (m^-1 sr^-1)
beta_R = @(R) beta_a.* exp(-R/8000); % backscatter profile (m^-1 sr^-1)

% Constants and conversions
c = 3e8;                    % Speed of light (m/s)
M_H2O = 18.015;             % Molar mass of water (g/mol)
N_A = 6.022e23;             % Avogadro's number (molecules/mol)
sigma_on = sigma_on_cm2 * 1e-4;    
sigma_off = sigma_off_cm2 * 1e-4;  
%delta_sigma = sigma_on - sigma_off; % Differential Absorption 
GM_to_NUM = N_A / M_H2O;    % (g/m^3) -> (molecules/m^3)
K_on = sigma_on * GM_to_NUM;
K_off = sigma_off * GM_to_NUM;
K_delta = K_on - K_off;      % Differential Absorption Factor

% Range Vector
R = (R_start : dR : R_end)';
N_bins = length(R);


%% Lidar Forward Model 

% System Constant, note 1e30 was empirically selected to get a reasonable number of summed counts per/bin
C_sys = (P_laser * c * tau_pulse / 2) * 1e30; 

% Geometric Overlap Function 
O_R_model = 1 - exp(-(R / O_scale).^2);
O_R_model(R < R_start) = 0; % Ensure 0 until the start of the retrieval
O_R = O_R_model;

% Calculate the Total Extinction Coefficient (alpha) in m^-1
alpha_off = alpha_non_abs + K_off * wv_mass_true(R);
alpha_on = alpha_non_abs + K_on * wv_mass_true(R);

% Calculate Optical Depth (OD) by summing alpha * dR
OD_off = cumsum( alpha_off * dR );
OD_on = cumsum( alpha_on * dR );

% True number online and offline counts
N_on_true = C_sys * (1./R.^2) .* O_R .* beta_R(R) .* exp(-2 * OD_on);
N_off_true = C_sys * (1./R.^2) .* O_R .* beta_R(R) .* exp(-2 * OD_off);


%% Detector Contamination Model

N_AP_peak = C_sys * AP_fraction;
N_AP_decay = N_AP_peak * exp(-(R - R_start) / lambda_AP);
N_AP_total = N_AP_decay + N_dark;

% Contaminated signal profiles 
N_off_raw = N_off_true + N_AP_total;
N_on_raw = N_on_true + N_AP_total;
N_bg = N_dark + N_solar;  

%% Contaminated Retrieval using analytical error addition

% calculate rho_retrieved = rho_true + rho_error
wv_retrieved = NaN(N_bins, 1);

Delta_R = dR;

for i = 1 : N_bins - 1
    % Get the true water vapor density
    R_mid = R(i) + Delta_R / 2;
    wv_true_at_mid = wv_mass_true(R_mid); % in g/m^3

    % Use N_off_true at the start of the bin for the error denominator
    N_off_true_at_R = N_off_true(i);

    % Calculate the AP derivative (slope) at R(i)
    R_for_slope = R(i);
    AP_slope = N_AP_peak * (-1/lambda_AP) * exp(-(R_for_slope - R_start) / lambda_AP);

    % Calculate the Error Term
    % rho_error ≈ -1/(2*K_delta) * (d(N_AP)/dR) / N_off_true
    % The negative sign on AP_slope and the external negative sign ensures a positive rho_error.
    wv_error = (-1 / (2 * K_delta)) * (AP_slope / N_off_true_at_R);

    % Final Contaminated Retrieval 
    wv_retrieved(i) = wv_true_at_mid + wv_error;
end


%% Error Calculation and Plotting

% Preparation for Plotting (N-1 length vectors) 
R_plot = R(1:end-1) + dR/2;
rho_true_plot = wv_mass_true(R_plot);
rho_retrieved_plot = wv_retrieved(1:end-1);

% Calculate the Systematic Relative Error 
Relative_Error = 100 * (rho_retrieved_plot - rho_true_plot) ./ rho_true_plot;

% Clamp the error magnitude for clean plotting
max_err = 500;
Relative_Error(Relative_Error > max_err) = max_err;
Relative_Error(Relative_Error < -max_err) = -max_err;

% Calculate Water Vapor Optical Depth at R_end (4 km) ---
R_4km_idx = find(R >= 4000, 1, 'first');
WV_OD_4km = (OD_on(R_4km_idx) - OD_off(R_4km_idx))% / 2;
disp('should be close to 1 for optimal profile retrievals');
disp('-------------------------------------------------------------------');


%% Plotting Results

figure('Position', [100, 100, 1400, 450]);

% Raw Measured Signals 
subplot(1, 3, 1);
semilogy(R, N_off_raw, 'b-', 'DisplayName', 'N_{off} Raw'); hold on;
semilogy(R, N_on_raw, 'r-', 'DisplayName', 'N_{on} Raw');
semilogy(R, N_AP_total, 'k--', 'LineWidth', 2, 'DisplayName', 'Afterpulse');
%semilogy(R, N_bg * ones(size(R)), 'g:', 'LineWidth', 1, 'DisplayName', 'Constant BG');
title('Raw Measured Signals (Semilogy)');
xlabel('Range (m)');
ylabel('Counts/Bin');
legend('Location', 'southwest');
xlim([R_start, R_end]);
ylim([N_dark/2, max(N_off_raw)*1.1]);
grid on;

% True vs. Retrieved Water Vapor Density (check that these are the same)
subplot(1, 3, 2);
plot(R_plot, rho_true_plot, 'k-', 'LineWidth', 3, 'DisplayName', 'wv_{true}'); hold on;
plot(R_plot, rho_retrieved_plot, 'r--', 'LineWidth', 2, 'DisplayName', 'wv_{retrieved}');
title('Water Vapor Profile Retrieval');
xlabel('Range (m)');
ylabel('Density \rho (g/m^3)');
legend('Location', 'northeast');
ylim([0, wv_mass_surf * 2]);
xlim([R_start, R_end/2]);
grid on;

% Systematic Relative Error 
subplot(1, 3, 3);
plot(R_plot, Relative_Error, 'm-', 'LineWidth', 2); hold on;
plot([R_start, R_end], [0, 0], 'k--');
title('Systematic Relative Error');
xlabel('Range (m)');
ylabel('Relative Error (%)');
ylim([-50, 50]);
xlim([R_start, R_end/2]);
grid on;
