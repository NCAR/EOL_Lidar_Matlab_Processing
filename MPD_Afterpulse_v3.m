%% MicroPulse DIAL Afterpulse Error Simulation (Refactored)

clear; close all; clc;

%% Initialization and Unit Conversion

% System Parameters
tau_pulse = 200e-9;         % Pulse width (s)
R_start = 200;              % Start range for retrieval (m)
R_end = 5000;               % End range for retrieval (m)
dR = 7.5;                   % Range resolution (m/bin)
P_laser = 50e-3;            % Avg laser power (W) 
sigma_on_cm2 =0.99e-23;    % Online cross-section (cm^2/molecule)
sigma_off_cm2 = 7.0e-25;    % Offline cross-section (cm^2/molecule)

% Detector Afterpulse Parameters
lambda_AP_1 = 90;             % Afterpulse decay length 1 (m)
lambda_AP_2 = 2000;           % Afterpulse decay length 2 (m)
A_1 = 40;                     % Scaling factor for short decay component
A_2 = 1;                      % Scaling factor for long decay component
AP_fraction = 2e-16;         % Afterpulse peak as fraction of C_sys

% Geometric Overlap Parameters
R_full_overlap = 1000;       % Range at which O(R) approaches 1 (m)
%O_scale = 100;               % Scaling factor for the exponential rise 

% Simple Atmosphere Parameters 
wv_mass_surf = 20;        % Water Vapor Density Profile
wv_mass_true = @(R) wv_mass_surf .* exp(-R/2000); % Mass density profile (g/m^3)
alpha_non_abs = 1e-8;     % Total aerosol and molecular (non-absorbing) extinction (m^-1)
beta_a = 1e-7;            % general backscatter coefficient (m^-1 sr^-1)
beta_R = @(R) beta_a.* exp(-R/8000); % backscatter profile (m^-1 sr^-1)

% Constants and conversions
c = 3e8;                    % Speed of light (m/s)
M_H2O = 18.015;             % Molar mass of water (g/mol)
N_A = 6.022e23;             % Avogadro's number (molecules/mol)
sigma_on = sigma_on_cm2 * 1e-4;    
sigma_off = sigma_off_cm2 * 1e-4;  
GM_to_NUM = N_A / M_H2O;    % (g/m^3) -> (molecules/m^3)
K_on = sigma_on * GM_to_NUM;
K_off = sigma_off * GM_to_NUM;
K_delta = K_on - K_off;      % Differential Absorption Factor

% Range Vector
R = (R_start : dR : R_end)';
N_bins = length(R);

% Midpoint range vector for plotting (N_bins - 1 length)
R_plot = R(1:end-1) + dR/2;

%% Lidar Forward Model 

% System Constant
C_sys = (P_laser * c * tau_pulse / 2) * 1e17; 

% Geometric Overlap Function 
O_R_model = 1 - exp(-(R / R_full_overlap).^2); % 
O_R_model(R < R_start) = 0; % Ensure 0 until the start of the retrieval
O_R = O_R_model;

% Extinction and Optical Depth
alpha_off = alpha_non_abs + K_off * wv_mass_true(R);
alpha_on = alpha_non_abs + K_on * wv_mass_true(R);
OD_off = cumsum( alpha_off * dR );
OD_on = cumsum( alpha_on * dR );

% True number online and offline counts
N_on_true = C_sys * (1./R.^2) .* O_R .* beta_R(R) .* exp(-2 * OD_on);
N_off_true = C_sys * (1./R.^2) .* O_R .* beta_R(R) .* exp(-2 * OD_off);

% Calculate the derivative of the true signal for the analytical error model
dNoff_dR = gradient(N_off_true, dR);

%% Detector Contamination Model
N_AP_peak = C_sys * AP_fraction;

% Bi-exponential afterpulse function
N_AP_decay = N_AP_peak * (A_1*exp(-(R - R_start)/lambda_AP_1)+ (A_2*exp(-(R - R_start) / lambda_AP_2)));

% Contaminated signal profiles
N_off_raw = N_off_true + N_AP_decay;
N_on_raw = N_on_true + N_AP_decay;

%% 1. FULL Contaminated Retrieval (Most Accurate Result - No Approximation)

% a. Calculate the contaminated log ratio
Log_Ratio_Raw = log(N_off_raw ./ N_on_raw);

% b. Calculate the derivative (slope) using the finite difference method
dLogRatio_dR_Raw = gradient(Log_Ratio_Raw, dR);

% c. Apply the DIAL equation to get the full retrieval
wv_retrieved_full = (1 / (2 * K_delta)) * dLogRatio_dR_Raw;

% SET PRIMARY RETRIEVAL TO FULL RESULT
wv_retrieved = wv_retrieved_full;

%% 2. ANALYTICAL Error Calculation (Vectorized Full Quotient Rule)
% This recalculates the analytical error term (rho_error_analytical) for comparison

% Use the first N_bins - 1 elements for derivatives and plotting
R_comp = R(1:end-1);
N_off_true_at_R = N_off_true(1:end-1);
N_AP_at_R = N_AP_decay(1:end-1);
dNoff_dR_at_R = dNoff_dR(1:end-1);

% Analytical AP_slope (dN_AP/dR) vectorized
AP_slope_analytic = N_AP_peak * (A_1*(-1/lambda_AP_1) * exp(-(R_comp - R_start)/lambda_AP_1) + ...
                       (A_2*(-1/lambda_AP_2) * exp(-(R_comp - R_start) / lambda_AP_2)));

% Full Quotient Rule (Analytic Approximation of the Error)
% Numerator: (N_sig * dNAP/dR) - (NAP * dNsig/dR)
Numerator_Analytic = (N_off_true_at_R .* AP_slope_analytic) - (N_AP_at_R .* dNoff_dR_at_R);
% Denominator: (N_sig)^2
Denominator_Analytic = N_off_true_at_R.^2;

% rho_error_analytical ≈ -1/(2*K_delta) * d/dR [ N_AP/N_sig ]
rho_error_analytical = (-1 / (2 * K_delta)) * (Numerator_Analytic ./ Denominator_Analytic);


%% Error Comparison and Plotting Setup
% True Water Vapor at plotting points
rho_true_plot = wv_mass_true(R_plot);
rho_retrieved_plot = wv_retrieved(1:end-1);

% --- True Systematic Error (Result of Full Retrieval - True Profile) ---
rho_error_true = rho_retrieved_plot - rho_true_plot;

% --- Error of the Analytical Approximation 
%Taylor_Error_Difference = 100 * (rho_error_analytical - rho_error_true) ./ rho_error_true; 
% --- Error of the Analytical Approximation Relative to the True Signal ---
% Formula: 100 * (Difference between approximated error and true error) / (True Signal)
Taylor_Error_Difference = 100 * (rho_error_analytical - rho_error_true) ./ rho_true_plot;



% Calculate and Clamp Total Relative Error (based on the accurate wv_retrieved)
Relative_Error = 100 * (rho_retrieved_plot - rho_true_plot) ./ rho_true_plot;
max_err = 500;
Relative_Error(Relative_Error > max_err) = max_err;
Relative_Error(Relative_Error < -max_err) = -max_err;

% Calculate Water Vapor Optical Depth
R_4km_idx = find(R >= 4000, 1, 'first');
WV_OD_4km = (OD_on(R_4km_idx) - OD_off(R_4km_idx)) % one-way WV OD
disp('should be close to 1 for optimal profile retrievals');
disp('-------------------------------------------------------------------');

% Display the Taylor approximation error using R_plot (R_mid)
R_mid = R_plot; 
fprintf('\n--------------------------------------------------------------\n');
fprintf('MAX ERROR OF ANALYTICAL APPROXIMATION (NEAR R<500m): %.2f%%\n', max(abs(Taylor_Error_Difference(R_mid < 500))));
fprintf('MAX ERROR OF ANALYTICAL APPROXIMATION (FAR R>3000m): %.2f%%\n', max(abs(Taylor_Error_Difference(R_mid > 3000))));
fprintf('--------------------------------------------------------------\n');


%% Plotting Results (Use the accurate wv_retrieved)
figure('Position', [100, 100, 1400, 450]);
% Raw Measured Signals 
subplot(1, 3, 1);
semilogy(R, N_off_raw, 'b-', 'DisplayName', 'N_{off} Raw'); hold on;
semilogy(R, N_on_raw, 'r-', 'DisplayName', 'N_{on} Raw');
semilogy(R, N_AP_decay, 'k--', 'LineWidth', 2, 'DisplayName', 'Afterpulse');
title('Raw Measured Signals (Semilogy)');
xlabel('Range (m)');
ylabel('Counts/Bin');
legend('Location', 'southwest');
xlim([R_start, R_end]);
ylim([1, max(N_off_raw)*1.1]);
grid on;
% True vs. Retrieved Water Vapor Density
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