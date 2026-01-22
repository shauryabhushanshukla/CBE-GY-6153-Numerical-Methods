% Clear environment
clc
clearvars
disp('Environment cleared and ready for new calculations.');

% Simple calculations (scalars)
y = 5;
x = 10;

total = x + y;
mul = x * y;         % scalar multiplication
div = y / x;

fprintf('Sum: %d\n', total);
fprintf('Multiplication: %d\n', mul);
fprintf('Division: %f\n', div);

% Vectors and matrices
x_row = 1:10;        % 1 through 10 (row vector)
whos
x_col = x_row.';     % transpose -> column vector

x_lin = linspace(1,10,100);  % 100 evenly spaced points between 1 and 10

y_vec = [12 50 -8 -100 -200];

A = [1,2; 3,4; 5,6; 7,8];     % 4x2 matrix

A_plus2 = A + 2;
A_times2 = A * 2;

At = A.';                     % transpose (2x4)
AAt = A * At;                 % (4x2)*(2x4) -> 4x4

w = linspace(0,100,10);       % start, stop, number of points
w_square = w .^ 2;            % element-wise square

x2 = 22:100;                  % another 1D array
Matrices = [1,2,3; 4,5,6];    % 2x3 matrix

Aones = ones(3,1);            % 3x1 column of ones
B10 = zeros(10,1);            % 10x1 zeros (use zeros(10) for 10x10)
C = zeros(2,8);               % 2x8 zeros matrix
D = eye(3);                   % 3x3 identity
