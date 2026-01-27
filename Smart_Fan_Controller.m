%% Question
% Revision Challenge: The "Smart Fan" Controller
%
% Imagine you are monitoring a CPU's temperature over 60 seconds.
%
% The goal is to design a smart fan controller that adjusts fan speed
% based on the CPU temperature.
%
% PARAMETERS:
% Time (t): 0 to 60 seconds
% Sampling interval: 1 second
%
% Temperature Formula:
% Temp(t) = 25 + 50 * sin(t/20)^2
%
% CONTROL LOGIC:
% Rule A (Idle):
% If Temp < 40, Fan Speed = 20%
%
% Rule B (Active):
% If 40 <= Temp <= 70, Fan Speed = 1.5 * (Temp - 40) + 20

% Rule C (Emergency):
% If Temp > 70, Fan Speed = 100%
%
% OBJECTIVE:
% Create a function called fan_controller that takes the temperature
% vector as input and returns the corresponding fan speed vector.
%
% Plot the temperature and fan speed on the same figure using subplot.
%
% Calculate the average fan speed over the entire 60-second duration.


%% Smart Fan Controller Simulation
clc;            % Clear command window
clearvars;      % Clear workspace variables
close all;      % Close open figures

% 1. PARAMETERS
t = 0:1:60;     % Time from 0 to 60 seconds (1-second intervals)

% 2. TEMPERATURE FORMULA 
% We use the dot operator ( .^ ) to ensure every element in t is squared
temp = 25 + 50 * sin(t/20).^2;

% 3. FUNCTION CALL
% We "pass" our temp vector into the function and receive 'fan_speed' back
speed = fan_controller(temp);

% 4. OBJECTIVE CALCULATIONS
avg_speed = mean(speed);
fprintf('Average Fan Speed over 60s: %.2f%%\n', avg_speed);

% 5. VISUALIZATION (Subplots)
% subplot(rows, columns, plot_number)
figure('Name', 'CPU Thermal Management');

subplot(2,1,1); 
plot(t, temp, 'r', 'LineWidth', 1.5);
grid on;
ylabel('Temp (°C)');
title('CPU Temperature Profile');

subplot(2,1,2);
plot(t, speed, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Time (seconds)');
ylabel('Fan Speed (%)');
title('Automated Fan Response');


%% --- CUSTOM FUNCTION AREA ---
% Local functions must stay at the very bottom of the script.
function fan_speed = fan_controller(temp_in)
    % A. PRE-ALLOCATION
    % Create a vector of zeros the same size as input for efficiency
    fan_speed = zeros(size(temp_in));
    
    % B. LOOPING THROUGH THE DATA
    % We must check every single temperature reading one by one
    for i = 1:length(temp_in)
        
        current_t = temp_in(i); % Look at the "i-th" temperature
        
        % C. CONTROL LOGIC (The Ladder)
        if current_t < 40
            % Rule A: Idle
            fan_speed(i) = 20;
            
        elseif current_t >= 40 && current_t <= 70
            % Rule B: Active (Linear scaling)
            fan_speed(i) = 1.5 * (current_t - 40) + 20;
            
        else
            % Rule C: Emergency (Everything > 70)
            fan_speed(i) = 100;
        end
    end
end

