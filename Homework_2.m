%% HOMEWORK :#02
% Student Name: [Shaurya Shukla]
% Date: 2026-02-02

clear; clc; close all;

%% ========================================================================
%% PROBLEM 1: CLEAN METHANE COMBUSTION
%% ========================================================================

%% P1 Part (iii): Base Case Simulation
% Objective: Solve for flow rates given F1 = 1000 and F3 = 950.

% 1. Define Parameters
F1_in = 1000;          % Fresh Feed (mol/s)
F3_in = 950;           % Recycle Controller Setting (mol/s)
r2_in = 0.1;           % Ratio CH4/CO2 in Stream 2
r4_in = 0.001;         % Ratio CH4/CO2 in Stream 4

% 2. Define Matrix A and Vector b
% Order of x: [mCH4(2), mCO2(2), mCH4(3), mCO2(3), mCH4(4), mCO2(4), Rcomb]

A = zeros(7,7);
b = zeros(7,1);

% Eq 1: Combustor CH4 Balance
A(1,1) = -1; A(1,3) = 1; A(1,7) = -1;  b(1) = -F1_in;
% Eq 2: Combustor CO2 Balance
A(2,2) = -1; A(2,4) = 1; A(2,7) = 1;   b(2) = 0;
% Eq 3: Separator CH4 Balance
A(3,1) = 1; A(3,3) = -1; A(3,5) = -1;  b(3) = 0;
% Eq 4: Separator CO2 Balance
A(4,2) = 1; A(4,4) = -1; A(4,6) = -1;  b(4) = 0;
% Eq 5: Ratio Stream 2
A(5,1) = 1; A(5,2) = -r2_in;           b(5) = 0;
% Eq 6: Ratio Stream 4
A(6,5) = 1; A(6,6) = -r4_in;           b(6) = 0;
% Eq 7: Controller Stream 3 (F3 Fixed)
A(7,3) = 1; A(7,4) = 1;                b(7) = F3_in;

% 3. Solve and Display
x_base = A \ b;

fprintf('--- Problem 1, Part (iii) Results ---\n');
disp('Matrix A:'); disp(A);
disp('Vector b:'); disp(b);

% Format Output Vector: [mCH4(1)...mCH4(4), mCO2(1)...mCO2(4), Rcomb]
output_vec = [F1_in, x_base(1), x_base(3), x_base(5), ...
              0,     x_base(2), x_base(4), x_base(6), ...
              x_base(7)];

disp('Final Output Vector:');
fprintf('mCH4(1-4): %8.2f %8.2f %8.2f %8.2f\n', output_vec(1:4));
fprintf('mCO2(1-4): %8.2f %8.2f %8.2f %8.2f\n', output_vec(5:8));
fprintf('Rcomb:     %8.2f\n', output_vec(9));
disp('Unique Solution? Yes, matrix is full rank.');


%% P1 Part (iv): Parametric Study (Conversion vs F3)
% Objective: Vary F3 from 100 to 50,000 and plot conversion.

F3_range = linspace(100, 50000, 100); 
conv_results = zeros(1, 100);

for i = 1:100
    b_loop = b; 
    b_loop(7) = F3_range(i); % Update Controller
    
    x_loop = A \ b_loop;
    
    mCH4_3 = x_loop(3); % Recycle CH4
    Rcomb  = x_loop(7);
    
    % Conversion = Reacted / (Fresh + Recycle)
    conv_results(i) = (Rcomb / (F1_in + mCH4_3)) * 100;
end

figure(1);
plot(F3_range, conv_results, 'LineWidth', 2);
grid on;
title('P1 Part (iv): Methane Conversion vs Recycle Rate');
xlabel('Recycle Flow Rate F3 (mol/s)');
ylabel('Single Pass Conversion (%)');


%% P1 Part (v): Controller Re-design (Stream 4)
% Objective: Check if system is solvable when F4 is fixed instead of F3.

fprintf('\n--- Problem 1, Part (v) Results ---\n');

% Modify Matrix A
A_new = A;
A_new(7,:) = 0;         % Clear old controller
A_new(7,5) = 1;         % Set coeff for mCH4(4)
A_new(7,6) = 1;         % Set coeff for mCO2(4)

% Modify Vector b (Target F4 = 1000)
b_new = b;
b_new(7) = 1000;

disp('New Matrix A (Last row changed):');
disp(A_new);

if rcond(A_new) < 1e-12
    disp('!!! RESULT: MATRIX IS SINGULAR !!!');
    disp('Unique Solution? NO.');
    disp('Reason: Global Mass Balance (In = Out) means F4 is already fixed by F1.');
    disp('Adding a controller for F4 creates a redundant equation.');
else
    disp(A_new \ b_new);
end


%% ========================================================================
%% PROBLEM 2: CO2 MEMBRANE SEPARATION
%% ========================================================================

%% P2 Part (i) & (ii): Single Case Plot
fprintf('\n--- Problem 2, Part (ii) Results ---\n');

L = 100e-6;          % Length (m)
phi = 1e-7;          % Permeability
k_prime = 1e2;       % Absorption Rate
N = 100;             % Grid Points

[x_mem, p_mem] = solve_membrane(L, phi, k_prime, N);

figure(2);
plot(x_mem, p_mem, 'LineWidth', 2, 'Color', 'r');
grid on;
title(['P2 Part (ii): Pressure Profile (k'' = ', num2str(k_prime), ')']);
xlabel('Position x (m)'); ylabel('Partial Pressure (bar)');


%% P2 Part (iii): Surface Plot (Varying k')
fprintf('Generating Surface Plot for Part (iii)...\n');

k_vals = logspace(1, 3, 100);
p_surface = zeros(100, N);

for i = 1:100
    [x_mem, p_temp] = solve_membrane(L, phi, k_vals(i), N);
    p_surface(i, :) = p_temp;
end

figure(3);
[X_grid, K_grid] = meshgrid(x_mem, k_vals);
surf(X_grid, K_grid, p_surface);
shading interp; colorbar;
set(gca, 'YScale', 'log');
xlabel('Position (m)'); ylabel('k'' (log scale)'); zlabel('Pressure (bar)');
title('P2 Part (iii): Pressure vs Position and k''');
view(45, 30);


%% P2 Part (iv): Finding Minimum Thickness L
fprintf('\n--- Problem 2, Part (iv) Results ---\n');

k_target = 100;
tolerance = 1e-2;
L_guess = 10e-6;
step_size = 5e-6;
found = false;

for i = 1:1000
    [x_iv, p_iv] = solve_membrane(L_guess, phi, k_target, N);
    
    % Check gradient at outlet
    dx = x_iv(2) - x_iv(1);
    dp_dx = abs((p_iv(end) - p_iv(end-1)) / dx);
    
    if dp_dx < tolerance
        found = true;
        fprintf('SUCCESS: Minimum Thickness Found!\n');
        fprintf('L = %.2e meters (%.2f microns)\n', L_guess, L_guess*1e6);
        fprintf('Final Gradient: %.4f bar/m\n', dp_dx);
        break;
    end
    L_guess = L_guess + step_size;
end

if ~found
    disp('Could not converge on a thickness.');
end


%% LOCAL FUNCTIONS (Should be at the end of the script)

function [x, p] = solve_membrane(L, phi, k, N)
    x = linspace(0, L, N);
    h = x(2) - x(1);
    
    beta = 2 + (k * h^2) / phi;
    
    A = zeros(N,N);
    b = zeros(N,1);
    
    % Inlet BC (P = 9 bar)
    A(1,1) = 1; b(1) = 9;
    
    % Interior Points (Finite Difference)
    for i = 2:N-1
        A(i,i-1) = 1;
        A(i,i)   = -beta;
        A(i,i+1) = 1;
    end
    
    % Outlet BC (P = 0 bar) (Since x_CO2 = 0)
    A(N,N) = 1; b(N) = 0;
    
    p = A \ b;
end