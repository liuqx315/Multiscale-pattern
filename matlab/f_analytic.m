function dy = f_analytic(t, y)
% usage: dy = f_analytic(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% May 2012
% All Rights Reserved

% model parameters
global Pdata;
ep = Pdata.ep;
lam = -1/ep;

% form the ODE RHS
dy = lam*y + 1/(1+t^2) - lam*atan(t);

% end function
