function J = J_VanderPol(t, y)
% usage: J = J_VanderPol(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% model parameters
global Pdata;
ep = Pdata.ep;
ep2 = Pdata.ep2;

% extract variables
u = y(1);
v = y(2);

% form the ODE Jacobian
Juu = 0;
Juv = 1;
Jvu = -ep2/ep;
Jvv = (1-u^2)/ep;
J = [Juu, Juv; Jvu, Jvv];

% end function
