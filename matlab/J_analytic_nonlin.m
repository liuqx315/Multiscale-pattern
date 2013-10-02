function J = J_analytic_nonlin(t, y)
% usage: J = J_analytic_nonlin(t, y)
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------


% form the ODE Jacobian
%    f = (t+1)*exp(-y)  =>  J = -(t+1)*exp(-y)
J = -(t+1)*exp(-y);

% end function
