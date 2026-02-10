%% Lecture -2 (Homework-2) - Linear Algebra


function laminar_flow()
% Solves laminar flow between two parallel plates using finite differences

clc; close all;

%% Physical constants
mu   = 1e-3;     % viscosity [Pa.s]
V_up = 5e-5;     % upper plate velocity [m/s]
B    = 1e-3;     % plate separation [m]
dpdx = -1;       % pressure gradient [Pa/m]

%% Grid parameters
N  = 1000;           % number of interior grid points
dy = B/(N+1);        % grid spacing

% RHS constant from discretization
G = -(dy^2/mu)*dpdx;

%% Allocate matrix and RHS
A = spalloc(N, N, 3*N-2);   % tridiagonal sparse matrix
b = zeros(N,1);

%% First interior node (near bottom wall: v = 0)
A(1,1) = 2;
A(1,2) = -1;
b(1)   = G;

%% Interior nodes
for j = 2:N-1
    A(j,j-1) = -1;
    A(j,j)   = 2;
    A(j,j+1) = -1;
    b(j)     = G;
end

%% Last interior node (near moving top wall)
A(N,N-1) = -1;
A(N,N)   = 2;
b(N)     = G + V_up;

%% Solve system
v = A\b;

%% Grid for plotting
y = linspace(dy, B-dy, N);

%% Analytical solution
v_anal = (V_up/B).*y + (1/(2*mu))*dpdx.*(y.^2 - B*y);

%% Plot
figure;
plot(y, v_anal, 'LineWidth', 2); hold on;
plot(y, v, '.', 'MarkerSize', 6);
xlabel('y [m]');
ylabel('Velocity v [m/s]');
legend('Analytical', 'Finite Difference');
title('Laminar Couette–Poiseuille Flow');
grid on;

end
