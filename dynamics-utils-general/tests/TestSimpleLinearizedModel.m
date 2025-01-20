% Formulate a simple linearized model to test - Fully Linear System
clc; clear all; close all;
addpath("..");

A = rand(6);
B = randn(6, 3);
f_nl = @(t, y, uk) A*y + B*uk;
dfdu = @(t, y, uk) B;
dfdx = @(t, y, uk) A;

lm = LinearizedModel(f_nl, dfdx, dfdu, 6, 3);


t_min = 0;
t_max = 0.9;
%% Validation
x0_ref = randn(6,1);
uk_ref = randn(3, 1);

[A_k, B_k, c_k, xkp1_ref] = integrate_linearized(x0_ref, uk_ref, t_min, t_max, lm);

x0_pert = x0_ref+randn(6, 1);
uk_pert = uk_ref+randn(3, 1);

[~, ~, ~, xpert_final] = integrate_linearized(x0_pert, uk_pert, t_min, t_max, lm);

error = xpert_final - (A_k*x0_pert+B_k*uk_pert+c_k)


%% Compare to known analytic results
Ak_error = A_k-expm(A*t_max)

fun = @(tau) expm(A*(t_max-tau))*B;
Bk_error = integral(fun, t_min, t_max, 'ArrayValued', true)-B_k