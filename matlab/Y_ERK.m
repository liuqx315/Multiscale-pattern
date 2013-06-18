function [y,yerr] = Y_ERK(z, Fdata)
% usage: [y,yerr] = Y_ERK(z, Fdata)
%
% Inputs:  z = current guesses for [z1, ..., zs]
%          Fdata = structure containing extra information for evaluating F.
% Outputs: y = glued-together time-evolved solution
%          yerr = estimated solution error vector
%
% This routine takes as input the intermediate-time states (z) for a
% multi-stage ERK method, and pieces them together to form the time-evolved
% solution y(t_{n+1}). 
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% extract ERK method information from Fdata
B = Fdata.B;
[Brows, Bcols] = size(B);
s = Bcols - 1;
c = B(1:s,1);
b = (B(s+1,2:s+1))';

% check to see if we have coefficients for embedding
if (Brows > Bcols)
   b2 = (B(s+2,2:s+1))';
else
   b2 = b;
end

% get some problem information
[zrows,zcols] = size(z);
nvar = zrows;
if (zcols ~= s)
   error('Y_ERK error: z has incorrect number of stages');
end

% call f at our guesses
f = zeros(nvar,s);
for is=1:s
   t = Fdata.t + Fdata.h*c(is);
   f(:,is) = feval(Fdata.fname, t, z(:,is));
end

% form the solution
%    ynew = yold + h*sum(b(j)*fj)
y = Fdata.yold + Fdata.h*f*b;

% form the error estimate
yerr = Fdata.h*f*(b-b2);

% end of function
