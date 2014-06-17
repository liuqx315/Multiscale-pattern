function dy = f_bruss(t, y)
% usage: dy = f_bruss(t, y)
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
a = Pdata.a; 
b = Pdata.b; 
ep = Pdata.ep;

% form the ODE RHS
dy = [a - (y(3)+1)*y(1) + y(1)*y(1)*y(2);
      y(3)*y(1) - y(1)*y(1)*y(2);
      (b-y(3))/ep - y(3)*y(1)];

% end function
