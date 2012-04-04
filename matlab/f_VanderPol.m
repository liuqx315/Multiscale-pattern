function dy = f_VanderPol(t, y)
% usage: dy = f_VanderPol(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% model parameters
global Pdata;
ep = Pdata.ep;
ep2 = Pdata.ep2;

% extract variables
u = y(1);
v = y(2);

% form the ODE RHS
du = v;
dv = (v - v*u^2 - ep2*u)/ep;
dy = [du; dv];

% end function
