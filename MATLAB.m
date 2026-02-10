%% Lecture -1 of CBE-GY-6153

%% Clear Environment
% Reset the command window and workspace, close figures.
clc                     % clear command window
clearvars               % remove variables from workspace
close all               % close all figure windows
disp('Environment cleared and ready for new calculations.'); % status

%% Simple Scalar Calculations
% Basic arithmetic with scalars
y = 5;                  % scalar y
x = 10;                 % scalar x

total = x + y;          % sum
mul = x * y;            % multiplication
divv = y / x;           % division (rename from 'div' to avoid function name)

fprintf('Sum: %d\n', total);                     % show sum
fprintf('Multiplication: %d\n', mul);            % show product
fprintf('Division: %f\n', divv);                 % show quotient

%% Vectors, Matrices and Linear Algebra
% Create common vector and matrix types for examples

x_row = 1:10;           % row vector: [1 2 ... 10]
x_col = x_row.';        % column vector: transpose of x_row

x_lin = linspace(1,10,100);  % 100 evenly spaced points from 1 to 10

y_vec = [12 50 -8 -100 -200];% example row vector with mixed signs

A = [1,2; 3,4; 5,6; 7,8];    % 4x2 matrix
A_plus2 = A + 2;             % add scalar 2 to every element
A_times2 = A * 2;            % scale matrix by 2

At = A.';                    % transpose of A (2x4)
AAt = A * At;                % A * A' -> 4x4 symmetric matrix

w = linspace(0,100,10);      % 10 points between 0 and 100
w_square = w .^ 2;           % element-wise square

x2 = 22:100;                 % another 1D integer vector
Matrices = [1,2,3; 4,5,6];   % 2x3 example matrix

Aones = ones(3,1);           % 3x1 column of ones
B10 = zeros(10,1);           % 10x1 column of zeros
C = zeros(2,8);              % 2x8 zero matrix
D = eye(3);                  % 3x3 identity matrix

x = 1:2:10;                  % vector from 1 to 10 with step 2 -> [1 3 5 7 9]

%% Indexing Examples
% Demonstrate how to extract and modify elements and subarrays

A = [5,3,4.2; 8,9,0];        % 2x3 matrix example
val_22 = A(2,2);             % element at row 2, column 2 -> 9
sum_2row = A(2,2) + A(2,1);  % add elements in row 2: 9 + 8

A = [1 2 3 4 5 6 7 8];       % 1x8 row vector
first = A(1);                % first element -> 1
second = A(2);               % second element -> 2
last = A(end);               % last element -> 8

A = [1 2 3; 4 5 6];          % 2x3 matrix
A(1,1) = 100;                % modify element (row1,col1) -> now 100
second_row = A(2,:);         % all columns of row 2 -> [4 5 6]
subset_2 = A(2,1:3);         % columns 1 to 3 of row 2 -> same as above
subset_2_to_end = A(2,1:end);% same as previous (useful with unknown size)

%% QUESTION 1: Function Analysis (0 < x < 5)
% Analyze y = -(x - 3).^2 + 10 on the interval (0,5)

X = linspace(0,5,100);           % dense sampling of x in [0,5]
Y = -(X - 3).^2 + 10;            % element-wise square with vector X

[Max_value, I] = max(Y);         % Max value and index in Y
X_max_value = X(I);              % x at which max occurs

Min_value = min(Y);              % minimum y on sampled grid
Imin = find(Y == Min_value, 1);  % index of first minimum
X_min_value = X(Imin);           % x at which min occurs

% Anonymous function for exact evaluation at any x
Yfun = @(x) -(x - 3).^2 + 10;    % use dot when x may be vectorized

% Evaluate function at an arbitrary x
y_at_20_2 = Yfun(20.2);          % demonstrates function outside sampled range

fprintf('Max Y = %.4f at X = %.4f\n', Max_value, X_max_value);
fprintf('Min Y = %.4f at X = %.4f\n', Min_value, X_min_value);
fprintf('Y(20.2) = %.4f\n', y_at_20_2);

%% SECTION 2: Curve Transformations and Plotting (-10 to 10)
% Compare vertical and horizontal shifts of the parabola

x = -10:10;                      % integer x from -10 to 10
y1 = -(x - 3).^2 + 10;           % original curve
y2 = -(x - 3).^2 + 15;           % vertical shift up by +5
y3 = -(x - 5).^2 + 10;           % horizontal shift right by +2

figure;                          % new figure window
plot(x, y1, 'm-*', 'LineWidth', 0.8); hold on
plot(x, y2, 'b--o', 'LineWidth', 0.8);
plot(x, y3, 'g:+', 'LineWidth', 0.8);
xlabel('X'); ylabel('Y');
title('Curve Transformations: vertical and horizontal shifts');
legend('y = -(x-3)^2 +10','y = +(15)','y = (x-5)','Location','best');
grid on;
hold off;

%% SECTION 3: Logic with sin(x)
% Find what percent of sin(x) values exceed 0.8 on x in [0,10]

x = linspace(0,10,1000);         % high-resolution x
y = sin(x);                      % sine values
y_check = 0.8;                   % threshold

figure;
plot(x,y,'-'); hold on
plot([0,10],[y_check,y_check],'r-','LineWidth',1); % horizontal threshold line
xlabel('x'); ylabel('sin(x)');
title('sin(x) and threshold y = 0.8');
grid on;
hold off;

y_greater = y > y_check;                         % logical array
final_percent = 100 * sum(y_greater) / length(y); % percent of points > 0.8

fprintf('Percent of y > %.2f is %.3f%%\n', y_check, final_percent);

% Use size instead of non-existent height/width functions
[rows, cols] = size(y);
fprintf('y vector dimensions: %d x %d\n', rows, cols);

%% SECTION 4: Random Numbers and Logic
% A) Generate 10 random integers from 1 to 5, count number of 3s.
% B) Display 'wow!' if at least 30% (i.e., 3 out of 10) are 3s.

A = randi(5,1,10);            % 1x10 vector of integers 1..5
fprintf('A = [%s]\n', sprintf('%d ', A));  % print the vector

count_threes = sum(A == 3);   % count how many elements equal 3
if count_threes >= 3          % check for at least 3 occurrences
    disp('wow!')              % display message
end

%% IF-ELSE Logic (User Input)
% Prompt user for a temperature and respond with a message

temp = input('Enter a temperature (numeric):\n');

if temp > 30
    disp('It''s a hot day.');        % note doubled single quote for apostrophe
elseif temp < 15
    disp('It''s a cold day.');
else
    disp('It''s a pleasant day.');
end

%% For Loop Example
% Demonstrates a fixed-count loop with formatted output
for k = 1:5
    fprintf('This is the repetition number: %d\n', k);
end

%% While Loop Example
% Countdown style loop using a condition
energy = 10;                   % initial battery percent
while energy > 0
    fprintf('Still Working....Battery Percent is %d\n', energy);
    energy = energy - 1;       % decrement battery each loop
end
disp('Out of energy');

%% Custom Local Function
% Local functions must appear at the end of a script file (supported in modern MATLAB)
my_answer = calc_power(4,5);   % call helper function
fprintf('4^5 = %d\n', my_answer);

function result = calc_power(base, power)
    % calc_power  Compute base^power with a warning for zero exponent
    %
    % Inputs:
    %   base  - numeric scalar base
    %   power - nonnegative integer exponent
    %
    % Output:
    %   result - base raised to power
    %
    if power == 0
        warning('Any number to the power 0 is 1.'); % small helpful message
    end
    result = base ^ power;
end


%% Practice Questions

% Energy Harvesting Simulation

% Imagine a sensor that generates power based on a sine wave
% (similar to a light sensor over a full day).
%
% The goal is to simulate the total energy collected over 24 hours.
%
% PARAMETERS:
% Time (t): 24 hours total
% Sampling interval: every 15 minutes (0.25 hours)
%
% Power Formula:
% P(t) = sin((pi * t) / 12) + 1
%
% SAFETY RULE:
% If the power P exceeds 1.8, the sensor enters "Overdrive" mode
% and only 80% of that power is stored.
%
% OBJECTIVE:
% Calculate the total energy collected by summing all power values
% over the 24-hour period (after applying the safety rule).

clc, clearvars, clear all
% Parameter
t =0:0.25:24;

% Power Formula 
P_raw = sin((pi * t)/ 12) + 1;

% Calling saftey function 
total_safe_power = apply_safety(P_raw);

% objective 
total_energy = sum(total_safe_power) * 0.25 %
fprintf("Total Energy collected: %d", total_energy)

% Plotting 
plot(t,P_raw, 'r--','LineWidth',2.5);
hold on;
plot(t,total_safe_power, 'b','LineWidth',3)
xlabel("Time")
ylabel("Power (w)")
title('Energy Harvesting with Safety Rule');
legend('Raw Power', 'Safe Power')
grid on;



% Function: Apply Safety Rule
function safe_power = apply_safety(raw_power)
    safe_power = raw_power; % % Make a copy of the input
    
% Apply the safety rule: if power > 1.8, store only 80% of it
    for i = 1:length(safe_power)
        if safe_power(i) > 1.8
            safe_power(i)  = safe_power(i) * 0.8;
        end
    end
end

