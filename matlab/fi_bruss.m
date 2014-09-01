function dy = fi_bruss(t, y)
% usage: dy = fi_bruss(t, y)
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
b = Pdata.b; 
ep = Pdata.ep;

% form the ODE RHS
dy = [0; 0; (b-y(3))/ep];

% end function
