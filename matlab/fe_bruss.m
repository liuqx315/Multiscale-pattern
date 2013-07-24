function dy = fe_bruss(t, y)
% usage: dy = fe_bruss(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% model parameters
global Pdata;
a = Pdata.a; 
b = Pdata.b; 

% form the ODE RHS
dy = [a - (y(3)+1)*y(1) + y(1)*y(1)*y(2);
      y(3)*y(1) - y(1)*y(1)*y(2);
     -y(3)*y(1)];

% end function
