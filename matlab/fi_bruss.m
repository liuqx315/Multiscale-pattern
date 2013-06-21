function dy = fi_bruss(t, y)
% usage: dy = fi_bruss(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% model parameters
global Pdata;
b = Pdata.b; 
ep = Pdata.ep;

% form the ODE RHS
dy = [0; 0; (b-y(3))/ep];

% end function
