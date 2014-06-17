function dt = EStab_bruss_1D(t, y)
% usage: dt = EStab_bruss_1D(t, y)
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% extract problem data
global Pdata;
d1 = Pdata.d1;
d2 = Pdata.d2;
dx = Pdata.dx;

% explicit stability based on diffusive CFL
dt = min(dx/dx/d1, dx/dx/d2);

% end function
