function dt = EStab_analytic(t, y)
% usage: dt = EStab_analytic(t, y)
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% model parameters
global Pdata;
ep = Pdata.ep;

% Jacobian of explicit part.  Since this is independent of u, just
% set to a huge number
dt = 1.0e30;

% end function
