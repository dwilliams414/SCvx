function [A_k, B_k, c_k, xkp1_ref] = cr3bp_discretization(xk, tk, tkp1, uk, mass_parameter)
%CR3BP_DISCRETIZATION: Obtain the A_k, B_k, c_k and xkp1_ref values for a
%trajectory under control in the CR3BP with Z.O.H.
arguments
    xk (6, 1) double
    tk (1, 1) double
    tkp1 (1, 1) double
    uk (3, 1) double
    mass_parameter (1, 1) double
end
    ode_fun = @(t, Y) cr3bp_discretized_dynamics(t, Y, uk, mass_parameter);

    odeopts = odeset('AbsTol', 1e-13, 'RelTol', 1e-13);

    stm0 = eye(6);
    Y0 = [xk; stm0(:); zeros(6*3+6, 1)];

    [t, x] = ode45(ode_fun, [tk tkp1], Y0, odeopts);

    xkp1_ref = x(end, 1:6)';
    A_k = reshape(x(end, 7:42), 6, 6);
    B_k = A_k * reshape(x(end, 42+1:42+6*3), 6, 3);
    c_k = A_k * x(end, end-5:end)';
end

