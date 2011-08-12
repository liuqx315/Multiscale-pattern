function J = Ji_VanderPol(t, y)
% usage: J = Ji_VanderPol(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% model parameters
global Pdata;
ep = Pdata.ep;

% extract variables
u = y(1);
v = y(2);

% form the ODE Jacobian
Juu = 0;
Juv = 0;
Jvu = -1/ep;
Jvv = (1-u^2)/ep;
J = [Juu, Juv; Jvu, Jvv];

% end function
