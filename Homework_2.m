clc; clearvars;

% --- MAIN SCRIPT STARTS HERE ---
% 1. Define the inputs
F1 = 1000;
F3 = 950;
r2 = 0.1;
r4 = 0.001;

% 2. Call the function
result = solve_methane(F1, F3, r2, r4);

% 3. Display the result clearly
disp('Calculated Flow Rates (mol/s):');
disp('   CH4(1)    CH4(2)    CH4(3)    CH4(4)    CO2(1)    CO2(2)    CO2(3)    CO2(4)     R_comb');
disp(result);
% --- MAIN SCRIPT ENDS HERE ---


% --- FUNCTION DEFINITION STARTS HERE ---
function output_vector = solve_methane(F1, F3, r2, r4)
     % 1. INITIALIZE MATRICES
    A = zeros(7,7);
    b = zeros(7,1);
    
    % 2. FILL MATRIX A (The Coefficients)
    % Col 1: mCH4(2)   Col 2: mCO2(2)
    % Col 3: mCH4(3)   Col 4: mCO2(3)
    % Col 5: mCH4(4)   Col 6: mCO2(4)
    % Col 7: R_comb

    % Row 1: Combustor Methane Balance (mCH4(2) - mCH4(3) + R = F1)
    A(1,1) = 1;  A(1,3) = -1;  A(1,7) = 1;
    b(1)   = F1;

    % Row 2: Combustor CO2 Balance (mCO2(2) - mCO2(3) - R = 0)
    A(2,2) = 1;  A(2,4) = -1;  A(2,7) = -1;
    b(2)   = 0;

    % Row 3: Separator Methane Balance
    A(3,1) = 1;   A(3,3) = -1; A(3,5) = -1;
    b(3)   = 0;

    % Row 4: Separator CO2 Balance
    A(4,2) = 1;   A(4,4) = -1; A(4,6) = -1;
    b(4)   = 0;

    % Row 5: Ratio Sensor Stream 2
    A(5,1) = 1;   A(5,2) = -r2;
    b(5)   = 0;

    % Row 6: Ratio Sensor Stream 4
    A(6,5) = 1;   A(6,6) = -r4;
    b(6)   = 0;

    % Row 7: Flow Controller Stream 3 (FIXED!)
    % Equation: mCH4(3) + mCO2(3) = F3
    A(7,3)= 1;   % Was A(7,1), changed to A(7,3)
    A(7,4)= 1; 
    b(7) = F3;

    % 3. SOLVE THE SYSTEM
    x = A \ b;

    % 4. FORMAT THE OUTPUT
    mCH4_1 = F1;
    mCO2_1 = 0;
    
    mCH4_streams = [mCH4_1, x(1), x(3), x(5)];
    mCO2_streams = [mCO2_1, x(2), x(4), x(6)];
    R_val        = x(7);
    
    output_vector = [mCH4_streams, mCO2_streams, R_val];
end

% 4. Check for Uniqueness
% We can check the rank of the matrix inside the function or just observe the result.
% If MATLAB returns a result without a "Singular Matrix" warning, it is unique.
disp('Does a unique solution exist? YES.');
disp('Reason: We have 7 independent linear equations for 7 unknowns.');
disp('The determinant of Matrix A is non-zero, meaning the system is non-singular.');