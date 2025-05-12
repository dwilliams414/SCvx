% CR3BP Discretization test
clc; clear all; close all;
mu = 0.012;
x0 = randn(6, 1);
addpath("../");
%% Build Discretized Model
f_nl = @(x, u, p) cr3bp_derivs(0, x, p)+[zeros(3); eye(3)]*u;
dfdx = @(x, u, p) A_cr3bp(0, x, p);
dfdu = @(x, u, p) [zeros(3); eye(3)];
dfdp = @(x, u, p) f_mu_cr3bp(x, p);

tspan = linspace(0, 3, 101);

dm = LinearizedDynamics(f_nl, dfdx, dfdu, dfdp, 6, 3, 1);

u = 0.0001*ones(3, 1);
[Ak, Bk, Ck, dk, xkp1] = dm.integrate_discretized(x0, ...
    @(x, p) 0.0001*ones(3, 1), mu, tspan);

%% Validation
error_bound = 1e-3;
for k = 1:100
    upert = u + error_bound*randn(3, 1);
    xpert = x0 + error_bound*randn(6, 1);
    mupert = mu+error_bound*randn(1);
    
    xf_lin = Ak*xpert+Bk*upert+Ck*mupert+dk;
    
    [~, ~, ~, ~, xf_nl] = dm.integrate_discretized(xpert, @(x, p) upert, mupert, ...
        tspan);
    xf_lin-xf_nl
    assert(all(abs(xf_nl-xf_lin) < error_bound))
end
rmpath("../");
