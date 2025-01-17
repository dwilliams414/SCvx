function [A_k, B_k, c_k, xkp1_ref] = integrate_linearized(xk, uk, tk, tkp1, lm, opts)
%INTEGRATE_LINEARIZED: Integrate linearized model to obtain discrete
%matrices
arguments
    xk (:, 1) double
    uk (:, 1) double
    tk (1, 1) double
    tkp1 (1, 1) double
    lm (1, 1) LinearizedModel
    opts.ode_solver (1, 1) function_handle = @ode113;
    opts.rel_tol (1, 1) double = 1e-12;
    opts.abs_tol (1, 1) double = 1e-12;
end
    dynamics = lm.dynamics4discretization();
    integration_fcn = @(t, Y) dynamics(t, Y, uk);

    stm0 = eye(lm.n_x);
    Y0 = [xk; stm0(:); zeros(lm.n_x*lm.n_u, 1); zeros(lm.n_x, 1)];


    intopts = odeset('AbsTol', opts.abs_tol, 'RelTol', opts.rel_tol);

    [t, Y] = opts.ode_solver(integration_fcn, [tk tkp1], Y0, intopts);

    A_k = reshape(Y(end, lm.stm_indices), [lm.n_x, lm.n_x]);
    B_k = A_k*reshape(Y(end, lm.B_indices), [lm.n_x lm.n_u]);
    c_k = A_k*reshape(Y(end, lm.c_indices), lm.n_x, 1);
    xkp1_ref = Y(end, 1:6)';
end

