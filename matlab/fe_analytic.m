function dy = fe_analytic(t, y)
% usage: dy = fe_analytic(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% form the ODE RHS
dy = 1/(1+t^2);

% end function
