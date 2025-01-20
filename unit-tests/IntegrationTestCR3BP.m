% Test Integration and Discretization
addpath("../dynamics-utils-cr3bp-mex/");
addpath("../dynamics-utils-general/");

%%  Create a Simple system to try out, CR3BP
clc; clear all; close all;
cr3bp_sys = load("em_constants.mat");
mass_parameter = cr3bp_sys.mu;
orb_dat = load("L2SHalo_FullStandardized.mat");
i_orb = 15;
x0 = orb_dat.x0_array(:, i_orb);
IP = orb_dat.IP_vals(i_orb)*0.8;
%% Specify Functions - CR3BP with No Control
uk = 1e-3*randn(3, 1);
f_nl = @(t, x, u) controlled_cr3bp(t, x, mass_parameter, @(t, x) uk);
dfdx = @(t, x, u) A_cr3bp(t, x, mass_parameter);
dfdu = @(t, x, u) [zeros(3); eye(3)];
n_x = 6;
n_u = 3;

lm = LinearizedModel(f_nl, dfdx, dfdu, n_x, n_u);

%% Specify Initial State and control

% Using Linearized Model implementation
xk = x0;
tic;
[A_k, B_k, c_k, xkp1_ref] = integrate_linearized(xk, uk, 0, IP, lm);
toc;

% Using Mex Implementation
tic;
[Akmex, Bkmex, ckmex, xkp1mex] = cr3bp_discretization_mex(xk, 0, IP, uk, mass_parameter);
toc;

tic;
[~, ~, ~, ~] = cr3bp_discretization(xk, 0, IP, uk, mass_parameter);
toc;

%% Errors in linearization implementation
Ak_error = max(abs(A_k-Akmex), [], 'all')
Bk_error = max(abs(B_k-Bkmex), [], 'all')
ck_error = max(abs(c_k-ckmex), [], 'all')
final_state_error_comparison = max(abs(xkp1mex-xkp1_ref), [], 'all')

%% Errors with respect to nonlinear
x0_pert = x0 + randn(6, 1)*1e-5;
uk_pert = uk + randn(3, 1)*1e-5;

ode_fun = @(t, y) controlled_cr3bp(t, y, mass_parameter, @(t, y) uk_pert);
odeopts = odeset('AbsTol', 1e-13, 'RelTol', 1e-13);
[tnl, xnl] = ode113(ode_fun, [0 IP], x0_pert, odeopts);

% Errors from Prediction
final_state_diff_actual = xnl(end, :)'-xkp1mex
final_state_diff_pred = Akmex*(x0_pert-x0)+Bkmex*(uk_pert-uk)

final_state_diff_error_nl = final_state_diff_actual-final_state_diff_pred

final_state_error_nl = xnl(end, :)'-(Akmex*x0_pert+Bkmex*uk_pert+ckmex)

%% Validate against STM - if control is set to zero
phi0 = eye(6);
[t, x] = IntegrateCR3BP([x0; phi0(:)], [0 IP], mass_parameter);
stm = state2phimats(x(end, :)');