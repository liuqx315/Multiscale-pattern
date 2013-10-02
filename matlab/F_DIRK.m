function F = F_DIRK(z, Fdata)
% usage: F = F_DIRK(z, Fdata)
%
% Inputs:  z = current guess for stage solution
%          Fdata = structure containing extra information for evaluating F.
% Outputs: F = residual at current guess
%
% This function computes the (non)linear residuals for an intermediate
% stage solution, through calling the user-supplied (in Fdata) ODE
% right-hand side function.
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% extract DIRK method information from Fdata
B = Fdata.B;
[Brows, Bcols] = size(B);
s = Bcols - 1;
c = B(1:s,1);
A = B(1:s,2:s+1);
h = Fdata.h;
st = Fdata.stage;
t = Fdata.t + Fdata.h*c(st);

% form the DIRK residual
%    F = z - rhs - h*(a(stage,stage)*fstage)
F = z - Fdata.rhs - h*A(st,st)*feval(Fdata.fname, t, z);

% end of function
