function J = J_analytic_nonlin(t, y)
% usage: J = J_analytic_nonlin(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved


% form the ODE Jacobian
%    f = (t+1)*exp(-y)  =>  J = -(t+1)*exp(-y)
J = -(t+1)*exp(-y);

% end function
