% =========================================================================
% Rb87 D2 LINE SPECTROSCOPY MODEL (HSRL & DAVLL)
% This script models the absorption profile of Rubidium 87 under 
% varying magnetic fields, from Zeeman to Paschen-Back regimes.
% =========================================================================
close all
clear all

% --- PHYSICAL CONSTANTS ---
I_val = 1.5;     % Nuclear spin of Rb87 (3/2)
J_g   = 0.5;     % Ground state electronic angular momentum (S 1/2)
J_e   = 1.5;     % Excited state electronic angular momentum (P 3/2)
muB   = 1.3996;  % Bohr Magneton in MHz/Gauss (Frequency shift per Gauss)

% --- HYPERFINE CONSTANTS (MHz) ---
% These represent the "coupling strength" between the nucleus and electron
A_g = 3417.3413; % Ground state magnetic dipole constant
gJ_g = 2.0023;   % Ground state Lande g-factor
A_e = 84.7185;   % Excited state magnetic dipole constant
gJ_e = 1.3362;   % Excited state Lande g-factor
gI = -0.000995;  % Nuclear g-factor

% --- CALIBRATION (Pinning the HSRL Lock Point) ---
% We define "0 MHz" as the strongest transition (F=2 -> F'=3) at low field.
% Eg(8) is the highest energy ground sublevel; Ee(16) is the highest excited.
[E0g, ~, ~, ~] = diagonalize_H(I_val, J_g, A_g, gJ_g, gI, 0.1, muB);
[E0e, ~, ~, ~] = diagonalize_H(I_val, J_e, A_e, gJ_e, gI, 0.1, muB);
f_ref = E0e(16) - E0g(8); 

% --- SIMULATION PARAMETERS ---
fields = [0, 100, 500, 1000, 5000];    % Magnetic field values in Gauss
detuning_axis = linspace(-20000, 20000, 10000); % Frequency range in MHz
doppler_width = 500; % MHz (Standard Doppler broadening at ~70 deg C)

figure('Color', 'w', 'Name', 'Rb87 Hyperfine Evolution');

for f_idx = 1:length(fields)
    B = fields(f_idx);
    
    % 1. BUILD & DIAGONALIZE HAMILTONIANS
    % This finds the specific energy levels (E) and wavefunctions (V) for this B-field
    [Eg, Vg, mJg_l, mIg_l] = diagonalize_H(I_val, J_g, A_g, gJ_g, gI, B, muB);
    [Ee, Ve, mJe_l, mIe_l] = diagonalize_H(I_val, J_e, A_e, gJ_e, gI, B, muB);
    
    spec_p = zeros(size(detuning_axis)); % Sigma+ spectrum (Blue shifted)
    spec_m = zeros(size(detuning_axis)); % Sigma- spectrum (Red shifted)
    
    % 2. COMPUTE TRANSITIONS
    % We loop through all 8 ground states (g) and 16 excited states (e)
    for g = 1:8
        for e = 1:16
            freq = (Ee(e) - Eg(g)) - f_ref; % Transition frequency relative to lock
            
            % 3. CALCULATE TRANSITION PROBABILITY (MATRIX ELEMENTS)
            % This is the "Overlap" logic. We check how much these states connect.
            s_p = 0; s_m = 0;
            for i = 1:8  % Component index for ground wavefunction
                for j = 1:16 % Component index for excited wavefunction
                    % Selection Rule: mI (nuclear spin) does not flip
                    if abs(mIg_l(i) - mIe_l(j)) < 0.1 
                        % Overlap (Weight) of the state components
                        weight = (Vg(i,g) * Ve(j,e))^2;
                        
                        % Polarization identification based on delta mJ
                        if abs(mJe_l(j) - mJg_l(i) - 1) < 0.1     % Delta mJ = +1
                            s_p = s_p + weight;
                        elseif abs(mJe_l(j) - mJg_l(i) + 1) < 0.1 % Delta mJ = -1
                            s_m = s_m + weight;
                        end
                    end
                end
            end
            
            % 4. GENERATE DOPPLER PROFILES
            % pop = 1/8 because atoms are distributed across ground states
            pop = 1/8; 
            term_p = pop * s_p * exp(-(detuning_axis-freq).^2/(2*doppler_width^2));
            term_m = pop * s_m * exp(-(detuning_axis-freq).^2/(2*doppler_width^2));
            
            spec_p = spec_p + term_p;
            spec_m = spec_m + term_m;
        end
    end
    
    % --- PLOTTING ---
    subplot(length(fields), 1, f_idx);
    plot(detuning_axis/1000, spec_p, 'b', 'LineWidth', 1.5); hold on;
    plot(detuning_axis/1000, spec_m, 'r', 'LineWidth', 1.5);
    grid on; xlim([-20 20]); ylim([0 0.75]);
    line([0 0], [0 0.5], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    title(sprintf('Magnetic Field = %.1f Gauss', B));
end
xlabel('Detuning from Lock Point (GHz)');

% =========================================================================
% DIAGONALIZATION FUNCTION
% Builds the Hamiltonian matrix using the coupled basis (Zeeman + Hyperfine)
% =========================================================================
function [E, V, mJ_l, mI_l] = diagonalize_H(I, J, A, gJ, gI, B, muB)
    dim = (2*I+1)*(2*J+1);
    H = zeros(dim);
    mJ_vals = -J:J; mI_vals = -I:I;
    [mJ_grid, mI_grid] = meshgrid(mJ_vals, mI_vals);
    mJ_l = mJ_grid(:); mI_l = mI_grid(:);
    
    for i = 1:dim
        % DIAGONAL ELEMENTS: Zeeman shift + Iz*Jz Hyperfine term
        H(i,i) = (gJ*mJ_l(i) + gI*mI_l(i))*muB*B + A*mI_l(i)*mJ_l(i);
        
        for j = 1:dim
            % OFF-DIAGONAL ELEMENTS: I+J- and I-J+ ladder operators
            % These terms "mix" the states to form the F manifolds
            if mI_l(i) == mI_l(j)+1 && mJ_l(i) == mJ_l(j)-1
                H(i,j) = (A/2) * sqrt(I*(I+1)-mI_l(j)*(mI_l(j)+1)) * ...
                                 sqrt(J*(J+1)-mJ_l(j)*(mJ_l(j)-1));
            elseif mI_l(i) == mI_l(j)-1 && mJ_l(i) == mJ_l(j)+1
                H(i,j) = (A/2) * sqrt(I*(I+1)-mI_l(j)*(mI_l(j)-1)) * ...
                                 sqrt(J*(J+1)-mJ_l(j)*(mJ_l(j)+1));
            end
        end
    end
    % Solve for eigenvalues (E) and eigenvectors (V)
    [V_raw, D] = eig(H);
    [E, idx] = sort(diag(D)); % Sort by energy to separate F=1 and F=2
    V = V_raw(:, idx);
end