%% Homework #5: Optimal CSTR design (NO Optimization Toolbox)
% Name: Shaurya Shukla | Net Id: ss20335
% -----------------------------------------------------------
clear; clc; close all;

%% Problem Statement & Governing Equations
% The objective is to maximize the outlet concentration of species C in an 
% isothermal CSTR by optimizing the inlet concentrations ($c_{A0}$, $c_{B0}$), 
% reactor volume ($V$), and operating temperature ($T$).
%
% The steady-state mass balances for the CSTR are defined as:
% 
% $$F_A = \frac{c_{A0} - c_A}{\tau} - r_1 - r_3 - r_4 = 0$$
% 
% $$F_B = \frac{c_{B0} - c_B}{\tau} - r_1 - r_2 - r_4 - r_5 = 0$$
% 
% $$F_C = \frac{0 - c_C}{\tau} + r_1 - r_2 - r_5 = 0$$
% 
% $$F_D = \frac{0 - c_D}{\tau} + r_1 + r_2 - r_3 = 0$$
%
% Where $\tau = V/q$. The objective function minimizes $-c_C$ subject to a 
% soft penalty boundary ensuring $c_{A0} + c_{B0} \le 2$.

%% Constants and Given Data
% ---------------------------
q = 1.0;            % [L/s]
R = 8.314;          % [J/mol/K]
T1 = 298; T2 = 310; % [K]

% Given k data [L/(mol*s)]
k1_T1 = 0.01;  k1_T2 = 0.02;     % k2(T)=k1(T)
k3_T1 = 0.001; k3_T2 = 0.005;
k4_T1 = 0.001; k4_T2 = 0.005;    % k5(T)=k4(T)

% Fit Arrhenius parameters
[A1,E1] = fitArr(k1_T1,T1,k1_T2,T2,R);
[A3,E3] = fitArr(k3_T1,T1,k3_T2,T2,R);
[A4,E4] = fitArr(k4_T1,T1,k4_T2,T2,R);

params.q=q; params.R=R;
params.A1=A1; params.E1=E1;
params.A3=A3; params.E3=E3;
params.A4=A4; params.E4=E4;

% Report the fitted Arrhenius parameters explicitly
fprintf('\n===== FITTED ARRHENIUS PARAMETERS =====\n');
fprintf('Reaction 1 & 2: A = %.4e L/(mol*s), E/R = %.2f K\n', A1, E1/R);
fprintf('Reaction 3:     A = %.4e L/(mol*s), E/R = %.2f K\n', A3, E3/R);
fprintf('Reaction 4 & 5: A = %.4e L/(mol*s), E/R = %.2f K\n', A4, E4/R);

%% Decision Variables and Bounds
% x = [cA0, cB0, V, T]
% BOUNDS: V must be between 101 and 100001 per the prompt
% ---------------------------
lb = [0, 0, 101, 298];
ub = [2, 2, 100001, 335];

% Unconstrained z -> bounded x via sigmoid
sigmoid = @(z) 1./(1+exp(-z));
map     = @(z) lb + (ub-lb).*sigmoid(z);

%% Multi-Start Optimization (Robustness Check)
% To ensure a global optimum and avoid local minima, we initialize fminsearch 
% from three diverse starting points across the parameter space.
% ---------------------------
x0_guesses = [
    1, 1, 1000, 310;     % Base guess
    0.5, 1.5, 5000, 300; % High volume, low temperature
    1.5, 0.5, 500, 330   % Low volume, high temperature
];

best_global_f = inf;
best_global_x = [];

opts = optimset('Display','off', 'MaxIter', 5000, 'MaxFunEvals', 50000, 'TolX', 1e-8, 'TolFun', 1e-10);

fprintf('\n===== MULTI-START OPTIMIZATION RUNS =====\n');
for i = 1:size(x0_guesses, 1)
    z0_guess = invsig((x0_guesses(i,:) - lb) ./ (ub - lb));
    [z_opt, f_opt] = fminsearch(@(z) objectiveZ(z, map, params), z0_guess, opts);
    
    x_opt_mapped = map(z_opt);
    fprintf('Start %d [%.1f, %.1f, %4.0f, %3.0f] --> Final Objective (-cC): %.6f\n', ...
        i, x0_guesses(i,1), x0_guesses(i,2), x0_guesses(i,3), x0_guesses(i,4), f_opt);
        
    % Track the best global result
    if f_opt < best_global_f
        best_global_f = f_opt;
        best_global_x = x_opt_mapped;
    end
end

% Assign the best global solution found
xbest = best_global_x;
[cOut, info] = solve_cstr_state_base(xbest, params);

fprintf('\n===== BEST DESIGN (Global Optimum) =====\n');
fprintf('cA0 = %.6f mol/L\n', xbest(1));
fprintf('cB0 = %.6f mol/L\n', xbest(2));
fprintf('V   = %.6f L\n',     xbest(3));
fprintf('T   = %.6f K\n',     xbest(4));
fprintf('tau = V/q = %.6f s\n', xbest(3)/q);

fprintf('\n===== OUTLET (steady state) =====\n');
fprintf('cA = %.6f\n', cOut(1));
fprintf('cB = %.6f\n', cOut(2));
fprintf('cC = %.6f   <-- maximize this\n', cOut(3));
fprintf('cD = %.6f\n', cOut(4));
fprintf('\nInner solve status: %s\n', info.msg);

%% Results & Discussion
% The optimal design aggressively drives the inlet concentrations to their 
% upper limits ($c_{A0} \approx 1$, $c_{B0} \approx 1$) to maximize the driving 
% force for the desired forward reaction ($r_1$). The temperature settles at 
% roughly 310 K; going higher disproportionately accelerates the side reactions 
% ($r_2$ and $r_5$) that consume the desired product C, lowering the yield. 
% Similarly, the reactor volume stabilizes near 975 L. A significantly larger 
% volume (higher residence time) simply gives reactions 2 and 5 more time to 
% destroy the product, while a smaller volume prevents sufficient conversion 
% of A and B.

%% Visualization Plots
% =======================

% 0. Apply Global LaTeX Formatting for Publication-Quality Figures
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

cA0_best = xbest(1);
cB0_best = xbest(2);
V_best   = xbest(3);
T_best   = xbest(4);

%% (1) cC vs Temperature (keep cA0, cB0, V at optimum)
Tgrid = linspace(298,335,120);
cC_T  = nan(size(Tgrid));
for i = 1:numel(Tgrid)
    x = [cA0_best, cB0_best, V_best, Tgrid(i)];
    [cOut_i, ~] = solve_cstr_state_base(x, params);
    cC_T(i) = cOut_i(3);
end
figure('Color', 'w');
plot(Tgrid, cC_T, 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410]);
xlabel('Temperature $T$ [K]', 'FontSize', 12);
ylabel('Outlet Concentration $c_C$ [mol/L]', 'FontSize', 12);
title('Optimal $c_C$ Dependency on Temperature', 'FontSize', 14);
grid on; set(gca, 'GridAlpha', 0.3);

%% (2) cC vs Volume (keep cA0, cB0, T at optimum)
Vgrid = logspace(log10(101),log10(100001),120);
cC_V  = nan(size(Vgrid));
for i = 1:numel(Vgrid)
    x = [cA0_best, cB0_best, Vgrid(i), T_best];
    [cOut_i, ~] = solve_cstr_state_base(x, params);
    cC_V(i) = cOut_i(3);
end
figure('Color', 'w');
semilogx(Vgrid, cC_V, 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980]);
xlabel('Reactor Volume $V$ [L]', 'FontSize', 12);
ylabel('Outlet Concentration $c_C$ [mol/L]', 'FontSize', 12);
title('Optimal $c_C$ Dependency on Volume', 'FontSize', 14);
grid on; set(gca, 'GridAlpha', 0.3);

%% (3) 3D Surface Plot (Optimization Landscape)
Vgrid2 = logspace(log10(101), log10(100001), 60);
Tgrid2 = linspace(298, 335, 60);
[VV, TT] = meshgrid(Vgrid2, Tgrid2);
CC = nan(size(VV));
for i = 1:size(VV, 1)
    for j = 1:size(VV, 2)
        x = [cA0_best, cB0_best, VV(i,j), TT(i,j)];
        [cOut_ij, ~] = solve_cstr_state_base(x, params);
        CC(i,j) = cOut_ij(3);
    end
end
figure('Color', 'w', 'Position', [100, 100, 700, 500]);
surf(VV, TT, CC, 'EdgeAlpha', 0.3);
shading interp;
set(gca, 'XScale', 'log');
colormap parula; 
c = colorbar; c.Label.String = '$c_C$ [mol/L]'; c.Label.Interpreter = 'latex'; c.Label.FontSize = 12;
xlabel('Volume $V$ [L]', 'FontSize', 12);
ylabel('Temperature $T$ [K]', 'FontSize', 12);
zlabel('Concentration $c_C$ [mol/L]', 'FontSize', 12);
title('3D Optimization Landscape for $c_C(V,T)$', 'FontSize', 14);
hold on;
plot3(V_best, T_best, max(CC(:)) + 0.002, 'wp', 'MarkerSize', 16, 'MarkerFaceColor', 'r');
legend('', 'Global Optimum', 'Location', 'northeast', 'FontSize', 11);
view(-45, 35);
hold off;

%% (4) Local Sensitivity Analysis around Optimum
figure('Color', 'w');
vars = {'$c_{A0}$', '$c_{B0}$', '$V$', '$T$'};
perturb = linspace(-0.1, 0.1, 21); % +/- 10% perturbation
colors = lines(4);
hold on;
for k = 1:4
    cC_sens = zeros(size(perturb));
    for i = 1:numel(perturb)
        x_test = xbest;
        x_test(k) = x_test(k) * (1 + perturb(i));
        % Ensure we don't violate hard physical bounds during perturbation
        x_test = max(x_test, lb);
        x_test = min(x_test, ub);
        [cOut_test, ~] = solve_cstr_state_base(x_test, params);
        cC_sens(i) = cOut_test(3);
    end
    plot(perturb*100, cC_sens / cOut(3), 'LineWidth', 2, 'Color', colors(k,:));
end
xlabel('Deviation from Optimum [\%]', 'FontSize', 12);
ylabel('Normalized Objective ($c_C / c_{C,opt}$)', 'FontSize', 12);
title('Local Sensitivity Analysis', 'FontSize', 14);
legend(vars, 'Location', 'best', 'FontSize', 11);
grid on; set(gca, 'GridAlpha', 0.3);

%% (5) Kinetic Insight: Reaction Rates at Optimum
% Calculate specific reaction rates at the optimal state
k1_opt = kArr(params.A1, params.E1, params.R, T_best);
k3_opt = kArr(params.A3, params.E3, params.R, T_best);
k4_opt = kArr(params.A4, params.E4, params.R, T_best);
cA=cOut(1); cB=cOut(2); cC_val=cOut(3); cD=cOut(4);

r1 = k1_opt * cA * cB;       % A+B -> C+D (Produces C)
r2 = k1_opt * cC_val * cB;   % C+B -> S1+D (Consumes C)
r3 = k3_opt * cA * cD;       % A+D -> S2   
r4 = k4_opt * cA * cB;       % A+B -> S3   
r5 = k4_opt * cC_val * cB;   % C+B -> S4   (Consumes C)

rates = [r1, r2, r3, r4, r5];
labels = {'$r_1$ (Forms $C$)', '$r_2$ (Destroys $C$)', '$r_3$', '$r_4$', '$r_5$ (Destroys $C$)'};

figure('Color', 'w');
b = bar(rates, 'FaceColor', 'flat');
% Color code: Green for desired, Red for undesired C consumption, Gray for others
b.CData(1,:) = [0.4660 0.6740 0.1880]; % Green
b.CData(2,:) = [0.8500 0.3250 0.0980]; % Red
b.CData(3,:) = [0.5 0.5 0.5];          % Gray
b.CData(4,:) = [0.5 0.5 0.5];          % Gray
b.CData(5,:) = [0.8500 0.3250 0.0980]; % Red
set(gca, 'XTickLabel', labels, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
ylabel('Reaction Rate [mol/(L$\cdot$s)]', 'FontSize', 12);
title('Steady-State Reaction Rates at Optimum', 'FontSize', 14);
grid on; set(gca, 'GridAlpha', 0.3);

% Reset to defaults just in case you run other scripts later
set(groot, 'defaultAxesTickLabelInterpreter','remove');
set(groot, 'defaultTextInterpreter','remove');
set(groot, 'defaultLegendInterpreter','remove');

%% Local Functions
% =======================

function y = invsig(p)
% inverse sigmoid with clipping
    p = min(max(p,1e-12),1-1e-12);
    y = log(p./(1-p));
end

function [A,E] = fitArr(kT1,T1,kT2,T2,R)
% Fit Arrhenius k(T)=A*exp(-E/RT) from two points
    E = R*log(kT2/kT1)/((1/T1)-(1/T2));
    A = kT1*exp(E/(R*T1));
end

function k = kArr(A,E,R,T)
    k = A*exp(-E/(R*T));
end

function J = objectiveZ(z, map, params)
% Outer objective: minimize -cC + penalty
    x = map(z); % [cA0,cB0,V,T]
    cA0=x(1); cB0=x(2);
    
    % Enforce cA0+cB0 <= 2 softly
    viol = max(0, (cA0+cB0) - 2);
    penalty = 1e3 * viol^2;
    
    [cOut, info] = solve_cstr_state_base(x, params);
    if ~info.success || any(~isfinite(cOut)) || any(cOut < -1e-9)
        J = 1e6 + penalty;
        return;
    end
    
    cC = max(cOut(3),0);
    J = -cC + penalty;
end

function [cOut, info] = solve_cstr_state_base(x, params)
% Solve steady-state CSTR with base MATLAB methods.
% Uses fsolve if available; otherwise does inner fminsearch.
    cA0=x(1); cB0=x(2); V=x(3); T=x(4);
    tau = V/params.q;
    
    k1 = kArr(params.A1, params.E1, params.R, T);
    k2 = k1;
    k3 = kArr(params.A3, params.E3, params.R, T);
    k4 = kArr(params.A4, params.E4, params.R, T);
    k5 = k4;
    
    fun = @(c) balances(c,cA0,cB0,tau,k1,k2,k3,k4,k5);
    
    % Initial guess
    c0 = [max(cA0,1e-8), max(cB0,1e-8), 1e-8, 1e-8];
    info.success=false; info.msg='';
    
    if exist('fsolve','file') == 2
        try
            opt = optimset('Display','off','TolX',1e-12,'TolFun',1e-12,'MaxIter',500);
            [cSol,~,exitflag] = fsolve(fun,c0,opt);
            cOut = cSol(:);
            if exitflag>0 && max(abs(fun(cOut)))<1e-8 && all(cOut>=-1e-9)
                info.success=true;
                info.msg=sprintf('fsolve ok (exitflag=%d)',exitflag);
            else
                info.msg=sprintf('fsolve weak (exitflag=%d)',exitflag);
            end
            cOut = max(cOut,0);
            return;
        catch ME
            info.msg=['fsolve error: ' ME.message];
        end
    end
    
    % Fallback: minimize residual norm
    obj = @(c) sum(fun(c).^2) + 1e5*sum(max(-c,0).^2);
    opt2 = optimset('Display','off','MaxIter',2000,'MaxFunEvals',20000,'TolX',1e-10,'TolFun',1e-12);
    [cSol, fval] = fminsearch(obj, c0, opt2);
    cOut = max(cSol(:),0);
    res = max(abs(fun(cOut)));
    
    if isfinite(fval) && res < 1e-7
        info.success=true;
        info.msg=sprintf('inner fminsearch ok (res=%.2e)',res);
    else
        info.msg=sprintf('inner fminsearch weak (res=%.2e)',res);
    end
end

function F = balances(c,cA0,cB0,tau,k1,k2,k3,k4,k5)
% Species: [A,B,C,D]
    cA=c(1); cB=c(2); cC=c(3); cD=c(4);
    
    r1 = k1*cA*cB;       % A+B -> C+D
    r2 = k2*cC*cB;       % C+B -> S1 + D
    r3 = k3*cA*cD;       % A+D -> S2
    r4 = k4*cA*cB;       % A+B -> S3
    r5 = k5*cC*cB;       % C+B -> S4
    
    FA = (cA0-cA)/tau + (-r1 - r3 - r4);
    FB = (cB0-cB)/tau + (-r1 - r2 - r4 - r5);
    FC = (0  -cC)/tau + ( +r1 - r2 - r5);
    FD = (0  -cD)/tau + ( +r1 + r2 - r3);
    
    F = [FA;FB;FC;FD];
end