function uxm = computeUxm(xvec, cr3bp_mu)
%computeUxm: Compute the partial derivative of Ux with respect to mass
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
    y = xvec(2);
    z = xvec(3);


    uxm = (x+2*cr3bp_mu-1)/d^3 + 3*(1-cr3bp_mu)*(x+cr3bp_mu)^2/d^5;
    uxm = uxm + (1-2*cr3bp_mu-x)/r^3 + 3*(cr3bp_mu*(x-1+cr3bp_mu)^2)/r^5;
end