%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cr3bp_derivs: Circular Restricted Three Body Problem Derivatives function
% Used for integration of trajectories in CR3BP
% Inputs:
% t: time vector
% state: state vector
%     x(1) = x (coordinate x)
%     x(2) = y
%     x(3) = z
%     x(4) = xdot (coordinate x)
%     x(5) = ydot
%     x(6) = zdot
% mu: CR3BP mu parameter
%
% Outputs:
% stateDeriv: derivative of state vector
% Validated 09/27/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function stateDeriv = cr3bp_derivs(t, state, mu)
   stateDeriv = zeros(6,1);
   
   % Define Useful Variables
   x = state(1);
   y = state(2);
   z = state(3);
   xdot = state(4);
   ydot = state(5);
   zdot = state(6);
   
   d = sqrt((x+mu)^2+y^2+z^2);
   r = sqrt((x-1+mu)^2+y^2+z^2);
   stateDeriv(1:3) = state(4:6);
   
   stateDeriv(4) = -(1-mu)*(x+mu)/d^3-mu*(x-1+mu)/r^3+x+2*ydot;
   stateDeriv(5) = -(1-mu)*y/d^3-mu*y/r^3+y-2*xdot;
   stateDeriv(6) = -(1-mu)*z/d^3-mu*z/r^3;
     
end

