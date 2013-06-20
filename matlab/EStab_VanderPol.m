function dt = EStab_VanderPol(t, y)
% usage: dt = EStab_VanderPol(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% model parameters
global Pdata;
ep = Pdata.ep;

% Jacobian of explicit part, J = [0, 1; 0 0], has all zero eigenvalues, so
% explicit method is unconditionally stable.  Set to a huge number
dt = 1/ep;

% end function
