clc, clearvars, clear all

A = [1 1 1; 2 1 3; 3 1 6]
b = [4; 7; 2]
X = A\b   % A * X  = b


% % Dot product of complex vector
% v1 = [1 + 2i; 2 - 1i; 3 + 0i]; % Define a complex vector
% v2 = [4 - 1i; 5 + 2i; 6 + 3i]; % Define another complex vector
% 
% result  = v1.^v2
% % dotProduct = dot(v1, v2); % Compute the dot product
% % disp(dotProduct)

v1 = [2, 3, 4];
v2 = [2, 3, 2];
result = v1.^v2 % Result: [4, 27, 16]

