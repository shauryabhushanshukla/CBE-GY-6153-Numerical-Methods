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

% Matrices, Arrays and Linear Algebra
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

x =1:2:10                     % start, jump, stop

% Index 
clc, clearvars

A = [5,3,4.2; 8 9 0]

A (2,2)    % for calling 9 means 2 row X 2 column

A(2,2) + A(2,1) % Adding 9 + 8

% If we are having matrix 

A = [1 2 3 4 5 6 7 8];

A(1)
A(2)

A(end) % Last value in the Matrix


A = [1 2 3; 4 5 6];
A(1,1) = 100

A(2,:)  % all the values of 2 row

A(2, 1:3)  % the values of 2 row from 1 column to 3 column 

A(2, 1:end)


%% QUESTION 1: Function Analysis
% A) What is the maximum value of the following equation on the range 0 < x < 5?
%    Equation: y = -(x - 3)^2 + 10
%
% B) What is the minimum of the function?
%
% C) At what x-value does the maximum y-value occur?
%
% D) What is y(20.7)?

clc, clearvars, close all  % close all clears all figures that we in upperscript


%generate bunch of X values 0 and 5

X = linspace(0,5)   % it will give 100 values

Y = (-(X-3).^2)  + 10   % its and element wise opeation add a dot 

%plot(X,Y, '*');

help max

[Max_value, I] = max(Y) % I will give the index of this 

X_max_value = X(I)   % 60th avalue of X


% Now we are going to define the custom fucntion, MAtlab calls these Anonymous Function 
Y = @(x) (-(x-3)).^2 +10
Y(20.2)

%% SECTION 2: Curve Transformations and Plotting (-10 to 10)

% A) Plot the equation y = -(x - 3)^2 + 10 from x = -10 to 10.

% B) How does the curve change if 15 is added instead of 10?

% C) How does the curve change if (x-5) is in parenthesis?
clc, clearvars

x = -10:10
y_1 = -(x-3).^2 +10
y_2 = -(x-3).^2 +15
y_3 = -(x-5).^2 +10

% figure(1) % Open one window
% plot(x, y_1, 'ms', 'LineWidth', 0.6)
% hold on    % Now everything after this stays in Figure 1
% plot(x, y_2, 'bv', 'LineWidth', 0.6)
% plot(x, y_3, 'g+', 'LineWidth', 0.6)
% 
% xlabel('X'), ylabel('Y') 
% legend('Y1', 'Y2', 'Y3')
% grid on


% xlim([0,2])    % Changing the limit of X from 0 to 2

%% Sublpots 


subplot(2,2,1)
plot(x,y_1,'*')
xlabel('X')
ylabel('Y')
title('Y vs X ')
legend('Original')

%% SECTION 3: LOGIC 
% A) Based upon the following equation: y = sin(x)
%    What percent of y-values are greater than 0.8 for x = 0 to 10?

clc, clearvars

% parameters
x = linspace(0,10,1000)
y =sin(x)
y_check = 0.8

% actions
plot(x,y,'*')

hold on 
plot([0,10],[0.8,0.8], 'r-')

y_greater  = y > y_check

Final_percent = sum(y_greater)/length(y)

height(y), width(y)





%%  LOGIC AND LOOPS
% IF approach 


% SECTION 4: Random Numbers and Logic
% A) Generate 10 random values from 1 to 5. Count the number of 3's.
% B) Display 'wow!' if more than 20% of the random numbers are 3.
% C) Do parts A and B with a For Loop.
% D) Extend to 10 million random numbers - which method is faster?



clc, clearvars

A = randi(5,1,10)
A
% e.g. A = [2 3 1 3 3 5 4 1 2 3]

if sum(A == 3) >= 3
    disp('wow!')
end


%% FOR LOOPS
    
clc, clearvars

for i  = 1:10
    i
end    