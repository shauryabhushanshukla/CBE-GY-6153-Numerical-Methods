% --- PART (i) ---

% 1. Define Parameters
N = 100;
W = 0.05;            % width [m]
dy = W / (N - 1);    % step size 
k = 0.6;             % Conductivity [W/mK]
I0 = 1000;           % Solar constant [W/m2]

% 2. Heat Generation 
% 50% absorption over width W
Q_val = (0.50 * I0) / W;   % Q = 10,000 [W/m3]

% 3. Initialize Matrix A and Vector b
A = zeros(N, N);
b = zeros(N, 1);

% 4. Build the Linear System
for j = 1:N
    if j == 1
        % Bottom wall - Adiabatic (-T1 + T2 = 0)
        A(j, j) = -1;
        A(j, j+1) = 1;
        b(j) = 0;        
        
    elseif j == N
        % Top Wall - Fixed Temp
        A(j, j) = 1;
        b(j) = 300;        
        
    else
        % Internal Nodes
        A(j, j-1) = -1;
        A(j, j) = 2;
        A(j, j+1) = -1;
        
        
        b(j) = (Q_val * dy^2) / k;  
    end
end

% 5. Solve the System
T = A \ b;  % Gaussian Elimination

% 6. Plot the result 
figure(1);
y = linspace(0, W, N);
plot(y, T, 'LineWidth', 2);
xlabel('Position Y [m]');      
ylabel('Temperature T [K]');
title('Part (i) Temperature Profile'); 
grid on;                        