function J = J_test(t, y)
% usage: J = J_test(t, y)
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% model parameters
I1 = 0.5;
I2 = 2;
I3 = 3;

% form the ODE RHS
J = [0, (I2-I3)/I1*y(3), (I2-I3)/I1*y(2);
     (I3-I1)/I2*y(3), 0, (I3-I1)/I2*y(1);
     (I1-I2)/I3*y(2), (I1-I2)/I3*y(1), 0];

% end function
