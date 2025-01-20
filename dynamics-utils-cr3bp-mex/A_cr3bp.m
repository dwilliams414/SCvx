%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A_cr3bp: Jacobian matrix for 3D CR3BP EOMs
% 
% Validated: 11/3/2021
% Inputs:
% t: time at which to evaluate
% state: 6-dimensional position-velocity state (6x1 vector)
% mu: Mass parameter
% 
% Outputs:
% A: Jacobian Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = A_cr3bp(t, state, mu)
    A11 = zeros(3);
    A12 = eye(3);
    
    % Compute UXX Terms
    pos = state(1:3);
    
    % Same partials
    uxx = computeUxx(pos, mu);
    uyy = computeUyy(pos, mu);
    uzz = computeUzz(pos, mu);
    
    % Mixed partials
    uxy = computeUxy(pos, mu);
    uxz = computeUxz(pos, mu);
    uyz = computeUyz(pos, mu);
    
    % The others
    uyx = uxy;
    uzx = uxz;
    uzy = uyz;
    
    UXX = [uxx uxy uxz;
           uyx uyy uyz;
           uzx uzy uzz];
       
    Omega = [0 2 0;
            -2 0 0;
             0 0 0];
         
    A = [A11 A12; UXX Omega];
    
end