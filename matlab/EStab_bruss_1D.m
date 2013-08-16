function dt = EStab_bruss_1D(t, y)
% usage: dt = EStab_bruss_1D(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% extract problem data
global Pdata;
d1 = Pdata.d1;
d2 = Pdata.d2;
dx = Pdata.dx;

% explicit stability based on diffusive CFL
dt = min(dx/dx/d1, dx/dx/d2);

% end function
