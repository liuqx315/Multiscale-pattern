function dt = EStab_CH_1D(t, y)
% usage: dt = EStab_CH_1D(t, y)
%
%   dy = -\partial_xx (c^2 \partial_xx y - y(y^2 - 1)), on [0,1]
%      u_x = 0     at x=0,x=1
%      u_xxx = 0   at x=0,x=1
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% extract model parameters, mesh size; set shortcut constants
global Pdata;
c = Pdata.c;
dx = Pdata.dx;
n = Pdata.n;
dx2 = dx*dx;
c2 = c*c/dx2/dx2;

% explicit portion handles anti-diffusion terms
%     - [y(i+1) - 2*y(i) + y(i-1)]/dx2
% so this should be unstable;  just set to a base diffusive CFL condition
dt = dx2;

% end function
