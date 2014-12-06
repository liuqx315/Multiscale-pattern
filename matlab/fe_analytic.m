function dy = fe_analytic(t, y)
% usage: dy = fe_analytic(t, y)
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% form the ODE RHS
dy = 1/(1+t^2);

% end function
