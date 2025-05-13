% Testing script for implementation of SCvx* Algorithm.  Goal is to hit a
% set of patch points states in the CR3BP, obtaining a minimum fuel
% trajectory.
clc; clear all; close all;

%% System Definition
em_sys = load("em_constants.mat");
orb_dat = load("L2SHalo_FullStandardized.mat");
x0 = orb_dat.x0_array(:, 12);
IP = orb_dat.IP_vals(12);
tspan = linspace(0, IP, 5);
[tref, xref] = IntegrateCR3BP(x0, tspan, em_sys.mu);

ref_states = xref(1:end-1, :)'+0.05*randn(6, 4);
ref_times = tref(1:end-1);
z0 = [ref_states(:); zeros(3*(length(ref_times)-1), 1)];
cvx_problem = CR3BPPosContinuityMinFuel(ref_states, ref_times, em_sys.mu);

%% Build Problem
problem = SCvx(cvx_problem);

sol = problem.solve(z0);