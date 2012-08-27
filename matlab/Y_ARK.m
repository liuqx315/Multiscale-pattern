function [y,yerr] = Y_ARK(z, Fdata)
% usage: [y,yerr] = Y_ARK(z, Fdata)
%
% Inputs:  z = current guesses for [z1, ..., zs]
%          Fdata = structure containing extra information for evaluating F.
% Outputs: y = glued-together time-evolved solution
%          yerr = estimated solution error vector
%
% This routine takes as input the intermediate-time states (z) for a
% multi-stage DIRK method, and pieces them together to form the time-evolved
% solution y(t_{n+1}). 
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% extract DIRK method information from Fdata
Bi = Fdata.B;
[Brows, Bcols] = size(Bi);
si = Bcols - 1;
ci = Bi(1:si,1);
bi = (Bi(si+1,2:si+1))';

% check to see if we have coefficients for embedding
if (Brows > Bcols)
   bi2 = (Bi(si+2,2:si+1))';
else
   bi2 = bi;
end


% extract ERK method information from Fdata
Be = Fdata.Be;
[Brows, Bcols] = size(Be);
se = Bcols - 1;
ce = Be(1:se,1);
be = (Be(se+1,2:se+1))';

% check to see if we have coefficients for embedding
if (Brows > Bcols)
   be2 = (Be(se+2,2:se+1))';
else
   be2 = be;
end

% ensure that methods have same number of internal stages
if (se ~= si)
   error('Y_ARK error: explicit and implicit method stage mis-match!');
end

% get some problem information
[zrows,zcols] = size(z);
nvar = zrows;
if (zcols ~= se)
   error('Y_ARK error: z has incorrect number of stages');
end

% call fe and fi at our guesses
fe = zeros(nvar,si);
fi = zeros(nvar,si);
for is=1:si
   t = Fdata.t + Fdata.h*ci(is);
   fi(:,is) = feval(Fdata.fname, t, z(:,is));
   t = Fdata.t + Fdata.h*ce(is);
   fe(:,is) = feval(Fdata.fnameE, t, z(:,is));
end

% form the solution
%    ynew = yol + h*sum_{j=1}^s (bi(j)*fi(zj) + be(j)*fe(zj))
y  = Fdata.yold + Fdata.h*fi*bi  + Fdata.h*fe*be;

% form the error estimate
yerr = Fdata.h*fi*(bi-bi2) + Fdata.h*fe*(be-be2);

% end of function
