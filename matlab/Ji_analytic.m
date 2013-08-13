function J = Ji_analytic(t, y)
% usage: J = Ji_analytic(t, y)
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

% form the ODE Jacobian
J = lam;

% end function
