function dy = f_analytic_nonlin(t, y)
% usage: dy = f_analytic_nonlin(t, y)
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% form the ODE RHS
%    f = (t+1)*exp(-y)
dy = (t+1)*exp(-y);

% end function
