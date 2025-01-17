% Create a Simple system to try out, prefereably one that we can validate
clc; clear all; close all;
cr3bp_sys = load("em_constants.mat");
mass_parameter = cr3bp_sys.mu;
orb_dat = load("L1L2Cycler_FullStandardized.mat");
x0 = orb_dat.x0_array(:, 301);

%% Specify Functions - CR3BP with No Control
f_nl = @(t, x, u) cr3bp_derivs(t, x, mass_parameter);
dfdx = @(t, x, u) A_cr3bp(t, x, mass_parameter);
dfdu = @(t, x, u) zeros(6, 3);
n_x = 6;
n_u = 3;

lm = LinearizedModel(f_nl, dfdx, dfdu, n_x, n_u);

%% Specify Initial State and control
xk = x0;
uk = zeros(3, 1);

[A_k, B_k, c_k, xkp1_ref] = integrate_linearized(xk, uk, 0, 0.2, lm);

% Compare with Nonlinear
xkp1 = A_k*xk+B_k*uk+c_k;

% Validate with Small perturbation
x0_new = xk + 1e-6*randn(6, 1);
[t_new, x_new] = IntegrateCR3BP(x0_new, [0 0.2], mass_parameter);
xkp1_new = x_new(end, :)';

xkp1_lin = A_k*x0_new+c_k