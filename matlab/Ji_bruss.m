function J = Ji_bruss(t, y)
% usage: J = Ji_bruss(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% model parameters
global Pdata;
ep = Pdata.ep;

% form the ODE Jacobian
J = [0, 0, 0; 0, 0, 0; 0, 0, -1/ep];

% end function
