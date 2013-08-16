function dt = EStab_analytic(t, y)
% usage: dt = EStab_analytic(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% May 2012
% All Rights Reserved

% model parameters
global Pdata;
ep = Pdata.ep;

% Jacobian of explicit part.  Since this is independent of u, just
% set to a huge number
dt = 1.0e30;

% end function
