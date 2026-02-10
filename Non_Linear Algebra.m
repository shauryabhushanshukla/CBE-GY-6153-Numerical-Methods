% --- PART (i): CONSTANT HEAT GENERATION ---
clear; clc;

% 1. Setup Parameters
N = 100;
W = 0.05;               % Width [m]
dy = W / (N - 1);       % Step size
k = 0.6;                % Conductivity
I0 = 1000;              % Solar Irradiance

% 2. Calculate Constant Q
Q_val = (0.50 * I0) / W; % Q = 10,000

% 3. Build Matrix System
A = zeros(N, N);
b = zeros(N, 1);

for j = 1:N
    if j == 1
        % Bottom: Adiabatic (-T1 + T2 = 0)
        A(j,j) = -1; A(j,j+1) = 1;
        b(j) = 0;
    elseif j == N
        % Top: Fixed Temp (TN = 300)
        A(j,j) = 1;
        b(j) = 300;
    else
        % Middle: Heat Equation
        A(j,j-1) = -1; A(j,j) = 2; A(j,j+1) = -1;
        b(j) = (Q_val * dy^2) / k;
    end
end

% 4. Solve and Plot
T_part1 = A \ b;

figure(1);
plot(linspace(0, W, N), T_part1, 'LineWidth', 2);
title('Part (i): Temperature Profile');
xlabel('Position y [m]'); ylabel('Temp [K]'); grid on;

% --- PART (ii): CALCULATE CONCENTRATION FACTOR ---

% 1. Get Average Temp from Part 1
T_avg_old = mean(T_part1);

% 2. Define Targets
T_target = 500;
T_ambient = 300;
C_old = 1;

% 3. Calculate Ratio
ratio = (T_target - T_ambient) / (T_avg_old - T_ambient);
C_new = ratio * C_old;

fprintf('Part (ii): Required Concentration Factor C = %.2f\n', C_new);




% --- PART (iii): NON-LINEAR ITERATIVE SOLVER ---
% clear; clc; % (Don't clear if you want to keep Part 1 results)

% 1. New Parameters
C = 300;                % Concentration Factor
I_in = C * 1000;        % Incoming Light (300,000)

% 2. Matrix A is the same as Part 1 (Geometry is constant)
% (We re-use Matrix A from Part 1 code above)

% 3. Initial Guess
T = ones(N, 1) * 300;   % Start with 300 K everywhere
error = 100;            % Big initial error
iter = 0;

% 4. The Loop
while error > 1e-4
    iter = iter + 1;
    T_old = T;
    
    % Step A: Update Epsilon
    epsilon = -1 - (1e-7 * T.^2);
    
    % Step B: Calculate Intensity I (Marching from Bottom)
    I = zeros(N, 1);
    I(1) = I_in;
    for j = 2:N
        % Formula: I_current = I_prev / (1 - dy * epsilon)
        I(j) = I(j-1) / (1 - dy * epsilon(j));
    end
    
    % Step C: Calculate Heat Q
    Q = -epsilon .* I;
    
    % Step D: Update Vector b
    b = zeros(N, 1);
    b(1) = 0;           % Bottom BC
    b(N) = 300;         % Top BC
    for j = 2:N-1
        b(j) = (Q(j) * dy^2) / k;
    end
    
    % Step E: Solve for New T
    T = A \ b;
    
    % Step F: Check Convergence
    error = max(abs(T - T_old));
end

fprintf('Part (iii): Converged in %d iterations.\n', iter);

% 5. Plot Results
figure(2);
subplot(2,1,1);
plot(linspace(0, W, N), T, 'r', 'LineWidth', 2);
title('Part (iii): Final Temperature Profile');
ylabel('Temp [K]'); grid on;

subplot(2,1,2);
plot(linspace(0, W, N), Q, 'b', 'LineWidth', 2);
title('Part (iii): Heat Generation Profile');
ylabel('Q [W/m^3]'); xlabel('Position y [m]'); grid on;


