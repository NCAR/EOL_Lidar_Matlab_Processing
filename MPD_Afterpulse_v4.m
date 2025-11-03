%% MicroPulse DIAL Afterpulse Error Simulation 

clear; close all; 

%% Initialization and Unit Conversion

% System Parameters
tau_pulse = 1e-6;           % Pulse width (s)
R_start = 200;              % Start range for retrieval (m)
R_end = 6000;               % End range for retrieval (m)
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
M = c*tau_pulse/2/dR; 

%% Lidar Forward Model 

% System Constant
C_sys = (P_laser * c * tau_pulse / 2) * C_sys_scale; 

% Geometric Overlap Function 
O_R_model = 1 ./ (1 + exp(-(R - R_full_overlap) / R_sigmoid_width));
O_R_model(R < R_start) = 0; % Ensure 0 until the start of the retrieval
O_R = O_R_model;
% figure(2)
%   semilogy(R, O_R)


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
N_AP_total = N_AP_decay + N_dark;

% Contaminated signal profiles
N_off_raw = N_off_true + N_AP_decay;
N_on_raw = N_on_true + N_AP_decay;

%% Full Numerical Retrieval of contaminated signals
% *** Applying 150 m (20 bin) smoothing for numerical stability ***



% Apply moving average filter to the raw signals
N_off_raw_smooth = movmean(N_off_raw, M);
N_on_raw_smooth = movmean(N_on_raw, M);

% Calculate the contaminated log ratio (using smoothed signals)
Log_Ratio_Raw = log(N_off_raw_smooth ./ N_on_raw_smooth);

% Calculate the derivative (slope) using the finite difference method
dLogRatio_dR_Raw = gradient(Log_Ratio_Raw, dR);

% Apply the DIAL equation to get the full retrieval
rho_retrieved_full = (1 / (2 * K_delta)) * dLogRatio_dR_Raw;

% Set Primary retrieval to the full result
rho_retrieved = rho_retrieved_full;

%% Analytical Error Calculation 
% This recalculates the analytical error term (rho_error_analytical) for comparison

% Term inside the derivative: N_AP * ( 1/N_off_true - 1/N_on_true )
Difference_Term = N_AP_decay .* (1 ./ N_off_true - 1 ./ N_on_true);

% Calculate the derivative of this term using gradient (numerical derivative)
dDifferenceTerm_dR = gradient(Difference_Term, dR);

% Apply the DIAL factor: -1/(2 * K_delta)
rho_error_analytical_full_R = (1 / (2 * K_delta)) * dDifferenceTerm_dR;


%% Error Comparison and Plotting Setup
% True Water Vapor at plotting points
rho_true_plot = wv_mass_true(R_plot);
rho_retrieved_plot = rho_retrieved(1:end-1);
rho_error_analytical = rho_error_analytical_full_R(1:end-1);

% True Systematic Error (Result of Full Retrieval - True Profile) 
rho_error_true = rho_retrieved_plot - rho_true_plot;

% Error of the Analytical Approximation Relative to the True Signal 
Analytical_Relative_Error = 100 * (rho_error_analytical) ./ rho_true_plot;

% % Calculate and Clamp Total Relative Error 
Relative_Error = 100 * (rho_retrieved_plot - rho_true_plot) ./ rho_true_plot;
% max_err = 500;
% Relative_Error(Relative_Error > max_err) = max_err;
% Relative_Error(Relative_Error < -max_err) = -max_err;

% Calculate Water Vapor Optical Depth
R_4km_idx = find(R >= 4000, 1, 'first');
WV_OD_4km = (OD_on(R_4km_idx) - OD_off(R_4km_idx)) % one-way WV OD
disp('should be close to 1 for optimal profile retrievals');
disp('-------------------------------------------------------------------');


%% Plotting Results 
figure('Position', [100, 100, 1400, 450]);

% Raw Measured Signals 
subplot(1, 4, 1);
semilogy(R, N_off_raw, 'k-', 'DisplayName', 'N_{off} Raw'); hold on;
semilogy(R, N_on_raw, 'b-', 'DisplayName', 'N_{on} Raw');
%semilogy(R, N_AP_total, 'g--', 'LineWidth', 2, 'DisplayName', 'Afterpulse+Dark');
semilogy(R, N_AP_decay, 'r--', 'LineWidth', 2, 'DisplayName', 'Afterpulse');
title('Raw Measured Signals');
xlabel('Range (m)');
ylabel('Counts/Bin');
legend('Location', 'southwest','FontSize', 10);
xlim([R_start, R_end]);
ylim([1, max(N_off_raw)*1.1]);
grid on;

% True vs. Retrieved Water Vapor Density
subplot(1, 4, 2);
plot(R_plot, rho_true_plot, 'k-', 'LineWidth', 3, 'DisplayName', 'wv_{true}'); hold on;
plot(R_plot, rho_retrieved_plot, 'b-', 'LineWidth', 2, 'DisplayName', 'wv_{numerical retrieval}');
title('Water Vapor Profile Retrieval');
xlabel('Range (m)');
ylabel('Density \rho (g/m^3)');
legend('Location', 'northeast','FontSize', 12);
ylim([0, wv_mass_surf * 2]);
xlim([0, R_end/2]);
grid on;

% Systematic density bias 
subplot(1, 4, 3);
plot(R_plot, rho_error_true,'b-', 'LineWidth', 2.5, 'DisplayName', 'Full Numerical' );
hold on
plot(R_plot, rho_error_analytical,'r--', 'LineWidth', 1.5, 'DisplayName', 'Analyitcal Approximation' );
%plot([R_start, R_end], [0, 0], 'k:', 'LineWidth', 1);
title('Systematic Bias');
xlabel('Range (m)');
ylabel('Systematic Density Bias (g/m^3)');
legend('Location','northeast','FontSize', 10);
ylim([-3, 3]);
xlim([0, R_end/2]);
grid on

% Systematic Relative Error 
 subplot(1, 4, 4);
 plot(R_plot, Relative_Error, 'b-', 'LineWidth', 2, 'DisplayName', 'Full Numerical' );
 hold on;
 plot(R_plot, Analytical_Relative_Error, 'r--', 'LineWidth', 2, 'DisplayName', 'Analyitcal Approximation' );
% plot([R_start, R_end], [0, 0], 'k--');
 title('Systematic Relative Error');
 xlabel('Range (m)');
 ylabel('Relative Error (%)');
 legend('Location','northeast','FontSize', 10);
 ylim([-20, 20]);
 xlim([0, R_end/2]);
 grid on;



figure('Position', [100, 100, 1000, 600]);

% Raw Measured Signals 
semilogy(R, N_off_raw, 'k-', 'DisplayName', 'N_{off} Raw'); hold on;
semilogy(R, N_on_raw, 'b-', 'DisplayName', 'N_{on} Raw');
%semilogy(R, N_AP_total, 'r--', 'LineWidth', 2, 'DisplayName', 'Afterpulse+Dark');
semilogy(R, N_AP_decay, 'r--', 'LineWidth', 2, 'DisplayName', 'Afterpulse');
title('Raw Measured Signals');
xlabel('Range (m)');
ylabel('Counts/Bin');
legend('Location', 'southwest','FontSize', 10);
xlim([0, R_end]);
ylim([1, max(N_off_raw)*1.1]);
grid on;