function J = Ji_bruss_1D(t, y)
% usage: J = Ji_bruss_1D(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% extract problem data
global Pdata;
a = Pdata.a; 
b = Pdata.b; 
d1 = Pdata.d1;
d2 = Pdata.d2;
m = Pdata.m;
dx = Pdata.dx;

% extract solution components
u = y(1:m);
v = y(m+1:2*m);

% initialize Jacobian blocks terms
Juu = sparse([],[],[],m,m,m);
Juv = sparse([],[],[],m,m,m);
Jvu = sparse([],[],[],m,m,m);
Jvv = sparse([],[],[],m,m,m);

% reaction components
for j=2:m-1
   Juu(j,j) = Juu(j,j) - (b+1) + 2*u(j)*v(j);
end
for j=2:m-1
   Juv(j,j) = Juv(j,j) + u(j)^2;
end
for j=2:m-1
   Jvv(j,j) = Jvv(j,j) - u(j)^2;
end
for j=2:m-1
   Jvu(j,j) = Jvu(j,j) + b - 2*u(j)*v(j);
end

% combine together into output
J = [Juu, Juv; Jvu, Jvv];

% end function
