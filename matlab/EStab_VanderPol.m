function dt = EStab_VanderPol(t, y)
% usage: dt = EStab_VanderPol(t, y)
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
ep = Pdata.ep;

% Jacobian of explicit part, J = [0, 1; 0 0], has all zero eigenvalues, so
% explicit method is unconditionally stable.  Set to a huge number
dt = 1/ep;

% end function
