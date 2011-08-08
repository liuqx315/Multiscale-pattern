function J = Jionization_2D(t, y, Pdata)
% usage: J = Jionization_2D(t, y, Pdata)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% model parameters
a = 0.6; 
b = 2; 
d1 = 0.025;
d2 = 0.025;

% determine spatial size, etc
n = length(y);
m = n/2;
T = y(1:m);
C = y(m+1:2*m);
dx = 1/(m-1);

% initialize Jacobian blocks terms
dTdT = sparse([],[],[],m,m,3*m);
dTdC = sparse([],[],[],m,m,m);
dCdT = sparse([],[],[],m,m,m);
dCdC = sparse([],[],[],m,m,3*m);

% diffusion components
for j=2:m-1
   dTdT(j,j-1) = dTdT(j,j-1) + d1/dx/dx;
   dTdT(j,j)   = dTdT(j,j) - 2*d1/dx/dx;
   dTdT(j,j+1) = dTdT(j,j+1) + d1/dx/dx;
end
for j=2:m-1
   dCdC(j,j-1) = dCdC(j,j-1) + d2/dx/dx;
   dCdC(j,j)   = dCdC(j,j) - 2*d2/dx/dx;
   dCdC(j,j+1) = dCdC(j,j+1) + d2/dx/dx;
end

% reaction components
for j=2:m-1
   dTdT(j,j) = dTdT(j,j) - (b+1) + 2*T(j)*C(j);
end
for j=2:m-1
   dTdC(j,j) = dTdC(j,j) + T(j)^2;
end
for j=2:m-1
   dCdC(j,j) = dCdC(j,j) - T(j)^2;
end
for j=2:m-1
   dCdT(j,j) = dCdT(j,j) + b - 2*T(j)*C(j);
end

% combine together into output
J = [dTdT, dTdC; dCdT, dCdC];

% end function
