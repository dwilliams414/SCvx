function dfdmu = f_mu_cr3bp(x, mu_cr3bp)
%f_mu_cr3bp: Compute the derivative of the system dynamics with respect to
%mass parameter
arguments
    x (6, 1) double
    mu_cr3bp (1, 1) double
end
    dfdmu = zeros(6, 1);
    dfdmu(4) = computeUxm(x, mu_cr3bp);
    dfdmu(5) = computeUym(x, mu_cr3bp);
    dfdmu(6) = computeUzm(x, mu_cr3bp);
end