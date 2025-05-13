% Testing script for implementation of SCvx* Algorithm.  Goal is to hit a
% set of patch points states in the CR3BP, obtaining a minimum fuel
% trajectory.
clc; clear all; close all;

%% System Definition
em_sys = load("em_constants.mat");

ref_times = [0 1 2];
ref_states = randn(6, 3);

cvx_problem = CR3BPPosContinuityMinFuel(ref_states, ref_times, em_sys.mu);