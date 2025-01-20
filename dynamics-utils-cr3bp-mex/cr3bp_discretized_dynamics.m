function dYdt = cr3bp_discretized_dynamics(t, Y, uk, mass_parameter)
%cr3bp_discretized_dynamics: CR3BP dynamics under acceleration control (uk)
%assuming Z.O.H.
arguments
    t (1, 1) double
    Y (:, 1) double
    uk (3, 1) double
    mass_parameter (1, 1) double
end
    n_x = 6;
    n_u = 3;

    % Indexing Variables
    state_indices = 1:n_x;
    stm_indices = state_indices(end)+1:state_indices(end)+n_x^2;
    ctrl_indices = stm_indices(end)+1:stm_indices(end)+n_x*n_u;
    c_vec_indices = ctrl_indices(end)+1:ctrl_indices(end)+n_x;

    dYdt = zeros(size(Y));
    
    % Obtain the B(t) matrix for control
    B = [zeros(3); eye(3)];

    % State Dynamics
    state = Y(state_indices);
    dYdt(state_indices) = cr3bp_derivs(t, state, mass_parameter) + ...
        B*uk;
    
    % STM Dynamics
    A_t = A_cr3bp(t, state, mass_parameter);
    stm = reshape(Y(stm_indices), 6, 6);
    dYdt(stm_indices) = reshape(A_t*stm, [], 1);

    % B_t Dynamics
    B_tau = B;
    dYdt(ctrl_indices) = reshape(inv(stm)*B_tau, [], 1);

    % c_k Dynamics
    c_t = cr3bp_derivs(t, state, mass_parameter) + B*uk - A_t*state - B_tau*uk;
    dYdt(c_vec_indices) = inv(stm)*c_t;
end