%% MicroPulse DIAL Afterpulse Error Simulation 

clear; close all; 

%% Initialization and Unit Conversion

% System Parameters
tau_pulse = 1e-6;           % Pulse width (s)
R_start = 200;                % Start range for retrieval (m)
R_end = 5000;               % End range for retrieval (m)
dR = 7.5;                   % Range resolution (m/bin)
P_laser = 50e-3;            % Avg laser power (W) 
sigma_on_cm2 =0.5e-23;      % Online cross-section (cm^2/molecule)
sigma_off_cm2 = 7.0e-25;    % Offline cross-section (cm^2/molecule)

% Detector Afterpulse Parameters
N_dark = 5;                   % Constant Dark Count Rate (counts/bin)
lambda_AP_1 = 80;            % Afterpulse decay length 1 (m)
lambda_AP_2 = 1000;           % Afterpulse decay length 2 (m)
A_1 = 1;                      % Scaling factor for short decay component
A_2 = 1/20;                   % Scaling factor for long decay component
AP_fraction = 4e-15;          % Afterpulse peak as fraction of C_sys
C_sys_scale = 1e17;         % Scales the counts/bin to something reasonable


% Geometric Overlap Parameters
R_full_overlap = 1200;       % Range at which O(R) approaches 1 (m)
R_sigmoid_width = 160;

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

% Smoothing to simulate laser pulse length
M = round(c*tau_pulse/2/dR); 
if M == 0; M = 1; end % Ensure M is at least 1 for movmean


%% Lidar Forward Model 

% System Constant
C_sys = (P_laser * c * tau_pulse / 2) * C_sys_scale; 

% Geometric Overlap Function 
O_R_model = 1 ./ (1 + exp(-(R - R_full_overlap) / R_sigmoid_width));
O_R_model(R < R_start) = 0; % Ensure 0 until the start of the retrieval
O_R = O_R_model;

% Extinction and Optical Depth
alpha_off = alpha_non_abs + K_off * wv_mass_true(R);
alpha_on = alpha_non_abs + K_on * wv_mass_true(R);
OD_off = cumsum( alpha_off * dR );
OD_on = cumsum( alpha_on * dR );

% True number online and offline counts (Instantaneous Model)
N_on_instantaneous = C_sys * (1./R.^2) .* O_R .* beta_R(R) .* exp(-2 * OD_on);
N_off_instantaneous = C_sys * (1./R.^2) .* O_R .* beta_R(R) .* exp(-2 * OD_off);

% Incorporate Pulse Smearing into the True Signal (N_on_true is M-bin smoothed)
N_on_true = movmean(N_on_instantaneous, M);
N_off_true = movmean(N_off_instantaneous, M);

%% Calculation of Baseline Retrieval (Error due to Pulse Smearing ONLY)
% N_off_true and N_on_true are the M-bin smoothed signals (N_lambda,baseline)
Log_Ratio_Baseline = log(N_off_true ./ N_on_true); 
dLogRatio_dR_Baseline = gradient(Log_Ratio_Baseline, dR);
rho_retrieved_baseline = (1 / (2 * K_delta)) * dLogRatio_dR_Baseline;
% The true water vapor profile (used for comparison)
rho_true_R_plot = wv_mass_true(R_plot);


%% Detector Contamination Model (Master Profile)

% This profile represents the maximum afterpulse measured during the 
% Zero-Signal Calibration (N_AP_zero).
N_AP_peak_zero = C_sys * AP_fraction;
N_AP_decay_zero = N_AP_peak_zero * (A_1*exp(-(R - R_start)/lambda_AP_1)+ (A_2*exp(-(R - R_start) / lambda_AP_2)));

% --- Define the True Residual Afterpulse in the Atmosphere (N_AP_true) ---
% Assume the true residual AP magnitude is a fixed fraction (e.g., 50%) of the max calibration profile.
AP_true_fraction = 0.5;
N_AP_true = AP_true_fraction * N_AP_decay_zero;

% Store the True Systematic Error for comparison (if no correction is applied)
N_off_raw_true = N_off_true + N_AP_true;
N_on_raw_true = N_on_true + N_AP_true;

% N_off/on_true is already M-bin smoothed. Use raw signal directly for log ratio.
Log_Ratio_True_Contaminated = log(N_off_raw_true ./ N_on_raw_true); 

rho_retrieved_true_error = (1 / (2 * K_delta)) * gradient(Log_Ratio_True_Contaminated, dR);

% Error is calculated relative to the smoothed baseline 
rho_error_true_uncorrected_AP = rho_retrieved_true_error(1:end-1) - rho_retrieved_baseline(1:end-1);


%% Scaling Sensitivity Test (Imperfect Subtraction)

% Scaling factor (kappa) now represents: N_AP_subtracted / N_AP_true
% A value of kappa=1.0 will yield the minimum retrieval error.
kappa_test_values = [0.25, 0.5, 0.75, 1.0, 1.25 1.5]; % Updated test values for clarity
error_profiles_scaled = cell(1, length(kappa_test_values));
legend_labels = cell(1, length(kappa_test_values));

for i = 1:length(kappa_test_values)
    kappa = kappa_test_values(i);
    legend_labels{i} = ['$\kappa = ' num2str(kappa) '$'];
    
    % Calculate the Afterpulse Profile to be SUBTRACTED (Scaled TRUE Profile)
    N_AP_sub = kappa * N_AP_true; 
    
    % Create the RAW Contaminated Signal (N_true + N_AP_true)
    N_off_raw = N_off_true + N_AP_true;
    N_on_raw = N_on_true + N_AP_true;
    
    % Create the CORRECTED Raw Signal (Contamination - Smoothed Subtraction)
    N_off_corrected = N_off_raw - N_AP_sub;
    N_on_corrected = N_on_raw - N_AP_sub;

    
    % Calculate the corrected log ratio (no extra smoothing needed here)
    Log_Ratio_Corrected = log(N_off_corrected ./ N_on_corrected);

    % Calculate the derivative and final retrieval
    dLogRatio_dR_Corrected = gradient(Log_Ratio_Corrected, dR);
    rho_retrieved_scaled = (1 / (2 * K_delta)) * dLogRatio_dR_Corrected;
    
    % Calculate the Systematic AP Error for this kappa relative to the smoothed baseline ---
    error_profiles_scaled{i} = rho_retrieved_scaled(1:end-1) - rho_retrieved_baseline(1:end-1);
end


%% Final Plotting of Sensitivity Results
rho_true_plot = wv_mass_true(R_plot);
disp('-------------------------------------------------------------------');
disp('True AP Fraction in Atmosphere (AP_true_fraction): '); disp(AP_true_fraction);
disp('Retrieval error will be minimized when kappa=1.0 (N_AP_subtracted / N_AP_true).');
disp('Error is now calculated relative to the pulse-smeared baseline retrieval.');
disp('Residual error at kappa=1.0 should be zero (or numerically negligible).');
disp('-------------------------------------------------------------------');


figure('Position', [100, 100, 500, 500]);

% Plot the corresponding Relative Error
% subplot(1, 2, 2);
hold on;
for i = 1:length(kappa_test_values)
    relative_error = 100 * error_profiles_scaled{i} ./ rho_true_plot;
    plot(R_plot, relative_error, 'LineWidth', 1.5, 'DisplayName', legend_labels{i});
end

% --- PLOT CHANGE ---
plot(R_plot, 100 * rho_error_true_uncorrected_AP ./ rho_true_plot, 'k--', 'LineWidth', 2.5, 'DisplayName', 'Uncorrected Baseline ($\kappa=0$ Subtraction)');

plot([R_start, R_end/2], [0, 0], 'k:');
title('Systematic Relative Error vs. Afterpulse Correction Factor ($\kappa$)', 'Interpreter', 'latex');
xlabel('Range (m)');
ylabel('Relative Error (%)');
legend('Location','northeast','FontSize', 10, 'Interpreter', 'latex');
ylim([-10, 10]);
xlim([R_start, R_end/2]);
grid on;


% figure('Position', [600, 100, 1000, 500]);
% 
% % Normalized O(R) on the left axis
% ax1 = subplot(1, 1, 1);
% yyaxis left;
% plot(R, O_R, 'b-', 'LineWidth', 2, 'DisplayName', 'Geometric Overlap O(R)');
% ylabel(ax1, 'Geometric Overlap O(R) (Normalized)');
% ylim([-0.1, 1.1]);
% xlim([R_start, R_end/2]);
% 
% % The K=1.0 Error Profile on the right axis
% yyaxis right;
% plot(R_plot, error_profiles_scaled{4}, 'r--', 'LineWidth', 2, 'DisplayName', 'Error at Kappa=1.0');
% ylabel(ax1, 'Systematic AP Error ($g/m^3$)');
% ylim([-0.05 0.05])
% 
% title('Systematic Afterpulse Error vs. Geometric Overlap Transition');
% xlabel('Range (m)');
% grid on;
% legend('Location','southeast','FontSize', 10);