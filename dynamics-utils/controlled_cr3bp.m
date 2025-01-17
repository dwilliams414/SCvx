function dxdt = controlled_cr3bp(t, x, mass_param, u)
%CONTROLLED_CR3BP CR3BP dynamics with control
    arguments
        t (1, 1) double
        x (6, 1) double
        mass_param (1, 1) double
        u (1, 1) function_handle = @(t, x) zeros(3, 1);  
    end

    dxdt = cr3bp_derivs(t, x, mass_param) + [zeros(3); eye(3)] * u(t, x);
end

