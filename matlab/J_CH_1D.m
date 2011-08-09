function J = J_CH_1D(t, y)
% usage: J = J_CH_1D(t, y)
%
%    J = -(c^2 \partial_{xxxx}(*) + (1 - 3y^2) \partial_{xx}(*))
%      u_x = 0     at x=0,x=1
%      u_xxx = 0   at x=0,x=1
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% extract model parameters, mesh size; set shortcut constants
global Pdata;
c = Pdata.c;
dx = Pdata.dx;
n = Pdata.n;
dx2 = dx*dx;
c2 = c*c/dx2/dx2;

% initialize Jacobian
J = sparse([],[],[],n,n,5*m);

% fourth-order terms
%   - c2*y(i+2)
%   + 4*c2*y(i+1) + y(i+1)^3/dx2 
%   - 6*c2*y(i) - 2*y(i)^3/dx2
%   + 4*c2*y(i-1) + y(i-1)^3/dx2
%   - c2*y(i-2)
for j=3:n-2
   J(j,j-2) = J(j,j-2) - c2;
   J(j,j-1) = J(j,j-1) + 4*c2 + 3/dx2*y(i-1)^2;
   J(j,j)   = J(j,j)   - 6*c2 - 6/dx2*y(i)^2;
   J(j,j+1) = J(j,j+1) + 4*c2 + 3/dx2*y(i+1)^2;
   J(j,j+2) = J(j,j+2) - c2;
end

% anti-diffusion terms
%     - [y(i+1) - 2*y(i) + y(i-1)]/dx2
for j=3:n-2
   J(j,j-1) = J(j,j-1) - 1/dx2;
   J(j,j)   = J(j,j)   + 2/dx2;
   J(j,j+1) = J(j,j+1) - 1/dx2;
end

% boundary terms
J(1,:) = 39/17*J(3,:) - 28/17*J(4,:) + 6/17*J(5,:);
J(2,:) = 67/34*J(3,:) - 21/17*J(4,:) + 9/34*J(5,:);
J(n-1,:) = 67/34*J(n-2,:) - 21/17*J(n-3,:) + 9/34*J(n-4,:);
J(n,:) = 39/17*J(n-2,:) - 28/17*J(n-3,:) + 6/17*J(n-4,:);


% end function
