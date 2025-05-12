function uzm = computeUzm(xvec, cr3bp_mu)
%computeUxm: Compute the partial derivative of Uz with respect to mass
%parmater.  Necessary for expanded STM that considers mu
% Validated 11/2/2023
arguments
    xvec (6, 1) double
    cr3bp_mu (1, 1) double
end
    d = compute_d(xvec(1:3), cr3bp_mu);
    r = compute_r(xvec(1:3), cr3bp_mu);

    % Unpack state
    x = xvec(1);
    % y = xvec(2);
    z = xvec(3);

    uzm = z/d^3+3*(1-cr3bp_mu)*z*(x+cr3bp_mu)/d^5;
    uzm = uzm - z/r^3 + 3*cr3bp_mu*z*(x-1+cr3bp_mu)/r^5;

end