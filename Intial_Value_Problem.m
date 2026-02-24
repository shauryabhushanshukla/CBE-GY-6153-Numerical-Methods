%% Numerical Methods in Chemical Engineering - HW#4
% *Name:* Shaurya Shukla | *Net Id:* ss20335
%
% This document models the spatiotemporal evolution of hydrogen and oxygen 
% concentration profiles within a Polymer Exchange Membrane (PEM) electrolyzer. 
%
% It investigates four operating scenarios:
%
% # Steady-state continuous operation.
% # Transient night shutdown with inert gas flushing.
% # Transient night shutdown with chamber mixing (no flushing).
% # Transient shutdown utilizing a reactive membrane to mitigate flammability.

clear; clc; close all;

%% 1. Physical Parameters and Initialization
% The geometric, thermodynamic, and transport properties are initialized based 
% on the system specifications. The system operates at 10 atm and 300 K. 
% Gases are assumed to behave ideally, and water content is neglected.

% Geometric Parameters
L      = 1e-4;           % Membrane thickness [m]
VoverA = 5e-3;           % Ratio of chamber volume to membrane area [m]

% Transport & Kinetic Parameters
DH2 = 1e-11;             % H2 diffusivity in membrane [m^2/s]
DO2 = 1e-12;             % O2 diffusivity in membrane [m^2/s]
k_L_per_mol_s = 1e-1;    % Rate constant for reactive membrane [L/mol-s]

% Thermodynamic States
P_atm = 10;              % Pressure [atm]
T = 300;                 % Temperature [K]
R = 8.314;               % Ideal gas constant [J/mol-K]
atm_to_Pa = 101325;      % Conversion factor [Pa/atm]

% Convert rate constant k to SI units: [m^3/mol-s]
k = k_L_per_mol_s * 1e-3;  

% Total gas concentration via Ideal Gas Law: C_tot = P/RT
P = P_atm * atm_to_Pa;   % Pressure [Pa]
Ctot = P/(R*T);          % Total concentration [mol/m^3]

% Discretization (Space)
N = 100;                 % Number of spatial nodes in membrane
x = linspace(0, L, N)';  % Column vector of positions
dx = x(2) - x(1);        % Spatial step size [m]

%% 2. Part (i): Steady-State Concentration Profiles
% Under continuous operation, boundary mole fractions are fixed. Assuming Fickian 
% diffusion at steady state, the profiles are strictly linear.

xH2_left  = 0.999;           % H2-evolution side
xH2_right = 0.01;            % O2 side (H2 crossover fraction)
xO2_left  = 1 - xH2_left;    % Assume only H2+O2
xO2_right = 1 - xH2_right;

% Boundary concentrations (Assuming partition coefficient = 1)
CH2_L = xH2_left  * Ctot;
CH2_R = xH2_right * Ctot;
CO2_L = xO2_left  * Ctot;
CO2_R = xO2_right * Ctot;

% Analytical linear steady-state profiles
CH2_ss = CH2_L + (CH2_R - CH2_L) * (x/L);  
CO2_ss = CO2_L + (CO2_R - CO2_L) * (x/L);

figure('Name','Part (i) Steady State Profiles');
plot(x, CH2_ss, 'b', 'LineWidth', 2); hold on;
plot(x, CO2_ss, 'r', 'LineWidth', 2);
xlabel('Position, x (m)'); ylabel('Concentration (mol/m^3)');
legend('H_2','O_2','Location','best');
title('Part (i): Steady-state concentration profiles (100 points)');
grid on;

%% 3. Numerical Formulation: Implicit Backward Euler Scheme
% An explicit Forward Time Centered Space (FTCS) scheme is subjected to strict 
% stability limits. To ensure absolute stability for extended simulations, 
% we implement an implicit Backward Euler finite difference scheme.
% 
% The discrete Laplacian operator (A) is formulated as a sparse tridiagonal matrix.

Ni = N - 2; % Number of interior nodes
e = ones(Ni,1);
A = spdiags([e -2*e e], -1:1, Ni, Ni);  % Tri-diagonal Laplacian stencil

% Helper anonymous function that advances one step for diffusion
diffuse_step_dirichlet = @(C, D, dt, BC_left, BC_right) ...
    stepDiffusionDirichlet(C, D, dt, dx, A, BC_left, BC_right);

%% 4. Part (ii): Night Shutdown with Inert Gas Flushing
% At t = 0 (night shutdown), the chambers are flushed with inert gas, 
% instantaneously dropping boundary concentrations to 0. We simulate the evolution 
% over the first 10 minutes.

t_end_ii = 10*60;              % 10 minutes = 600 s
dt_ii = 1.0;                   % 1 second time step
nt_ii = round(t_end_ii/dt_ii) + 1;
t_ii = linspace(0, t_end_ii, nt_ii);

CH2 = CH2_ss; CO2 = CO2_ss;
CH2_store = zeros(nt_ii, N); CO2_store = zeros(nt_ii, N);
CH2_store(1,:) = CH2'; CO2_store(1,:) = CO2';

BC0 = 0; % Boundary values during flushing are zero
for n = 2:nt_ii
    CH2 = diffuse_step_dirichlet(CH2, DH2, dt_ii, BC0, BC0);
    CO2 = diffuse_step_dirichlet(CO2, DO2, dt_ii, BC0, BC0);
    CH2_store(n,:) = CH2'; CO2_store(n,:) = CO2';
end

figure('Name','Part (ii) H2 Surface (Flushed)');
surf(x, t_ii, CH2_store, 'EdgeColor','none');
xlabel('x (m)'); ylabel('time (s)'); zlabel('C_{H2} (mol/m^3)');
title('Part (ii): H_2 concentration in membrane (flushed)');
view(35,30); colorbar; grid on;

figure('Name','Part (ii) O2 Surface (Flushed)');
surf(x, t_ii, CO2_store, 'EdgeColor','none');
xlabel('x (m)'); ylabel('time (s)'); zlabel('C_{O2} (mol/m^3)');
title('Part (ii): O_2 concentration in membrane (flushed)');
view(35,30); colorbar; grid on;

%% 5. Part (iii): Night Shutdown Without Flushing (Chambers Mix)
% Gases freely diffuse and mix in the chambers. The boundary concentrations 
% couple to the membrane flux via Fick's Law at the interfaces. We track the 
% system for 12 hours and monitor for flammability limits (0.04 < x_H2 < 0.96).

t_end_iii = 12*3600;           % 12 hours
dt_iii = 5.0;                  % 5 s timestep
nt_iii = round(t_end_iii/dt_iii) + 1;
t_iii = linspace(0, t_end_iii, nt_iii);

CH2 = CH2_ss; CO2 = CO2_ss;
CH2_Lch = CH2_L;   CO2_Lch = CO2_L;
CH2_Rch = CH2_R;   CO2_Rch = CO2_R;

stride = 10;  % Store every 10th timestep to optimize memory
idx_store = 1:stride:nt_iii;
t_store = t_iii(idx_store);
CH2_store_iii = zeros(length(idx_store), N);
CO2_store_iii = zeros(length(idx_store), N);

xH2_L = zeros(nt_iii,1); xH2_R = zeros(nt_iii,1);
flamm_low = 0.04; flamm_high = 0.96;
t_flamm_iii = NaN; store_k = 1;

for n = 1:nt_iii
    xH2_L(n) = CH2_Lch / (CH2_Lch + CO2_Lch);
    xH2_R(n) = CH2_Rch / (CH2_Rch + CO2_Rch);
    
    if isnan(t_flamm_iii) && ((xH2_L(n) > flamm_low && xH2_L(n) < flamm_high) || ...
                              (xH2_R(n) > flamm_low && xH2_R(n) < flamm_high))
        t_flamm_iii = t_iii(n);
    end
    
    if store_k <= length(idx_store) && n == idx_store(store_k)
        CH2_store_iii(store_k,:) = CH2'; CO2_store_iii(store_k,:) = CO2';
        store_k = store_k + 1;
    end
    if n == nt_iii, break; end

    % 1) Diffusion with Dirichlet BCs equal to chamber C
    CH2 = diffuse_step_dirichlet(CH2, DH2, dt_iii, CH2_Lch, CH2_Rch);
    CO2 = diffuse_step_dirichlet(CO2, DO2, dt_iii, CO2_Lch, CO2_Rch);
    
    % 2) Compute boundary fluxes
    JH2_left = +DH2 * (CH2(2) - CH2(1))/dx;
    JO2_left = +DO2 * (CO2(2) - CO2(1))/dx;
    JH2_right = -DH2 * (CH2(N) - CH2(N-1))/dx;
    JO2_right = -DO2 * (CO2(N) - CO2(N-1))/dx;
    
    % 3) Update well-mixed chambers
    CH2_Lch = max(CH2_Lch + (dt_iii/VoverA)*JH2_left, 0);
    CO2_Lch = max(CO2_Lch + (dt_iii/VoverA)*JO2_left, 0);
    CH2_Rch = max(CH2_Rch + (dt_iii/VoverA)*JH2_right, 0);
    CO2_Rch = max(CO2_Rch + (dt_iii/VoverA)*JO2_right, 0);
end

figure('Name','Part (iii) Chamber xH2 vs Time');
plot(t_iii/3600, xH2_L, 'b', 'LineWidth', 2); hold on;
plot(t_iii/3600, xH2_R, 'r', 'LineWidth', 2);
yline(flamm_low, '--k', '0.04'); yline(flamm_high, '--k', '0.96');
xlabel('Time (hours)'); ylabel('x_{H2} in chamber');
legend('Left chamber (H2 side)','Right chamber (O2 side)','Location','best');
title('Part (iii): Chamber H_2 mole fraction vs time'); grid on;

%% 6. Part (iv): Reactive Membrane Model
% We introduce a reaction term. We solve this using operator splitting: 
% an implicit diffusion step followed by local fixed-point iterations for 
% the non-linear reaction.

t_end_iv = 12*3600; dt_iv = 5.0;
nt_iv = round(t_end_iv/dt_iv) + 1;
t_iv = linspace(0, t_end_iv, nt_iv);

CH2 = CH2_ss; CO2 = CO2_ss;
CH2_Lch = CH2_L;   CO2_Lch = CO2_L;
CH2_Rch = CH2_R;   CO2_Rch = CO2_R;

t_store_iv = t_iv(1:10:nt_iv);
xH2_L_iv = zeros(nt_iv,1); xH2_R_iv = zeros(nt_iv,1);
t_flamm_iv = NaN;

for n = 1:nt_iv
    xH2_L_iv(n) = CH2_Lch / (CH2_Lch + CO2_Lch);
    xH2_R_iv(n) = CH2_Rch / (CH2_Rch + CO2_Rch);
    
    if isnan(t_flamm_iv) && ((xH2_L_iv(n) > flamm_low && xH2_L_iv(n) < flamm_high) || ...
                             (xH2_R_iv(n) > flamm_low && xH2_R_iv(n) < flamm_high))
        t_flamm_iv = t_iv(n);
    end
    if n == nt_iv, break; end
    
    % 1) Diffusion step (implicit)
    CH2 = diffuse_step_dirichlet(CH2, DH2, dt_iv, CH2_Lch, CH2_Rch);
    CO2 = diffuse_step_dirichlet(CO2, DO2, dt_iv, CO2_Lch, CO2_Rch);
    
    % 2) Reaction step (fixed-point iteration)
    CH2r = CH2; CO2r = CO2;
    for it = 1:6
        rate = k .* CH2r .* CO2r;   
        CH2r = max(CH2 - dt_iv * rate, 0);
        CO2r = max(CO2 - dt_iv * 0.5 * rate, 0);
    end
    CH2 = CH2r; CO2 = CO2r;
    
    % 3) Boundary fluxes and Chamber Update
    JH2_left  = +DH2 * (CH2(2) - CH2(1))/dx;
    JO2_left  = +DO2 * (CO2(2) - CO2(1))/dx;
    JH2_right = -DH2 * (CH2(N) - CH2(N-1))/dx;
    JO2_right = -DO2 * (CO2(N) - CO2(N-1))/dx;
    
    CH2_Lch = max(CH2_Lch + (dt_iv/VoverA)*JH2_left, 0);
    CO2_Lch = max(CO2_Lch + (dt_iv/VoverA)*JO2_left, 0);
    CH2_Rch = max(CH2_Rch + (dt_iv/VoverA)*JH2_right, 0);
    CO2_Rch = max(CO2_Rch + (dt_iv/VoverA)*JO2_right, 0);
end

figure('Name','Part (iv) Chamber xH2 vs Time (Reactive)');
plot(t_iv/3600, xH2_L_iv, 'b', 'LineWidth', 2); hold on;
plot(t_iv/3600, xH2_R_iv, 'r', 'LineWidth', 2);
yline(flamm_low, '--k', '0.04'); yline(flamm_high, '--k', '0.96');
xlabel('Time (hours)'); ylabel('x_{H2} in chamber');
legend('Left chamber (H2 side)','Right chamber (O2 side)','Location','best');
title('Part (iv): Chamber H_2 mole fraction vs time (Reactive)'); grid on;

% Console Output
fprintf('--- FLAMMABILITY SAFETY REPORT ---\n');
if isnan(t_flamm_iii)
    fprintf('Part (iii) Baseline: Safe. No crossover within 12 hrs.\n');
else
    fprintf('Part (iii) Baseline: Flammability limit reached at %.3f hrs.\n', t_flamm_iii/3600);
end
if isnan(t_flamm_iv)
    fprintf('Part (iv) Reactive: Safe. No crossover within 12 hrs.\n');
else
    fprintf('Part (iv) Reactive: Flammability limit reached at %.3f hrs.\n', t_flamm_iv/3600);
end

%% 7. Local Functions
function Cnew = stepDiffusionDirichlet(C, D, dt, dx, A, BC_left, BC_right)
    % Solves discrete implicit diffusion: (I - r*A) * C^{n+1} = C^n + BCs
    N = length(C);
    Ni = N - 2;                   
    r = D*dt/(dx^2);              
    
    M = speye(Ni) - r*A;
    b = C(2:N-1);
    
    b(1)   = b(1)   + r*BC_left;
    b(end) = b(end) + r*BC_right;
    
    Cnew = C;
    Cnew(1) = BC_left;
    Cnew(2:N-1) = M \ b;
    Cnew(N) = BC_right;
end