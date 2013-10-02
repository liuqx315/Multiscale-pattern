function dy = f_bruss_1D(t, y)
% usage: dy = f_bruss_1D(t, y)
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
a = Pdata.a; 
b = Pdata.b; 
d1 = Pdata.d1;
d2 = Pdata.d2;
m = Pdata.m;
dx = Pdata.dx;

% extract solution components
u = y(1:m);
v = y(m+1:2*m);

% initialize RHS terms
du = zeros(m,1);
dv = zeros(m,1);

% enforce stationary boundary conditions
du(1) = 0;  du(m) = 0;  dv(1) = 0;  dv(m) = 0;

% diffusion components
du(2:m-1) = d1/dx/dx*(u(3:m)+u(1:m-2)-2*u(2:m-1));
dv(2:m-1) = d2/dx/dx*(v(3:m)+v(1:m-2)-2*v(2:m-1));

% reaction components
du(2:m-1) = du(2:m-1) + a - (b+1)*u(2:m-1) + u(2:m-1).*u(2:m-1).*v(2:m-1);
dv(2:m-1) = dv(2:m-1) + b*u(2:m-1) - u(2:m-1).*u(2:m-1).*v(2:m-1);

% combine together into output
dy = [du; dv];

% end function
