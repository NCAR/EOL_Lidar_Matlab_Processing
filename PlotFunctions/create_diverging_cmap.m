function cmap = create_diverging_cmap(N, T_max, T_gray)
% CREATE_DIVERGING_CMAP Generates a custom diverging colormap (Blue-Gray-Red)
% with a fixed gray center that corresponds to a user-defined data range.
%
%   N:       Total number of colors in the map (e.g., 256).
%   T_max:   The maximum absolute data value (e.g., 5 for +/- 5K).
%   T_gray:  The absolute data value defining the extent of the gray center 
%            (e.g., 2 for +/- 2K).

% --- Define Key Colors ---
blue = [0 0 1];
red = [1 0 0];
light_gray = [0.9 0.9 0.9];

% --- Determine Indices for the Gray Band ---
% Calculate the ratio of the gray band range to the total range
ratio_gray = T_gray / T_max; 
if ratio_gray >= 1
    error('T_gray must be less than T_max.');
end

% N_gray is the total number of indices that will be dedicated to light_gray.
N_gray = round(N * ratio_gray);
if mod(N_gray, 2) ~= 0
    N_gray = N_gray + 1; % Ensure N_gray is an even number
end

N_half = N / 2;
N_gray_half = N_gray / 2;

% N_transition is the number of indices that will handle the smooth fade
N_transition = N_half - N_gray_half;

% --- Build the First Half (Blue -> Transition -> Light Gray) ---
% Part A: Blue to Gray Transition
R1_trans = linspace(blue(1), light_gray(1), N_transition)';
G1_trans = linspace(blue(2), light_gray(2), N_transition)';
B1_trans = linspace(blue(3), light_gray(3), N_transition)';
cmap_half1_trans = [R1_trans, G1_trans, B1_trans];

% Part B: Solid Gray Center (Lower Half)
cmap_half1_gray = repmat(light_gray, N_gray_half, 1);

% Concatenate the first half: Transition + Solid Gray
cmap_half1 = [cmap_half1_trans; cmap_half1_gray];

% --- Build the Second Half (Light Gray -> Transition -> Red) ---
% Part A: Solid Gray Center (Upper Half)
cmap_half2_gray = repmat(light_gray, N_gray_half, 1);

% Part B: Gray to Red Transition
R2_trans = linspace(light_gray(1), red(1), N_transition)';
G2_trans = linspace(light_gray(2), red(2), N_transition)';
B2_trans = linspace(light_gray(3), red(3), N_transition)';
cmap_half2_trans = [R2_trans, G2_trans, B2_trans];

% Concatenate the second half: Solid Gray + Transition
cmap_half2 = [cmap_half2_gray; cmap_half2_trans];

% --- Final Colormap ---
cmap = [cmap_half1; cmap_half2];

end
