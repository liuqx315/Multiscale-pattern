function dy = f_analytic_nonlin(t, y)
% usage: dy = f_analytic_nonlin(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

% form the ODE RHS
%    f = (t+1)*exp(-y)
dy = (t+1)*exp(-y);

% end function
