function dy = fCH_1D(t, y)
% usage: dy = fCH_1D(t, y)
%
%   dy = -\partial_xx (c^2 \partial_xx y - y(y^2 - 1)), on [0,1]
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

% initialize RHS
dy = zeros(n,1);

% fourth-order terms
%   - c2*y(i+2)
%   + 4*c2*y(i+1) + y(i+1)^3/dx2 
%   - 6*c2*y(i) - 2*y(i)^3/dx2
%   + 4*c2*y(i-1) + y(i-1)^3/dx2
%   - c2*y(i-2)
dy(3:n-2) = -c2*y(5:n) + 4*c2*y(4:n-1) + y(4:n-1).^3/dx2 ...
   - 6*c2*y(3:n-2) - 2*y(3:n-2).^3/dx2 + 4*c2*y(2:n-3) ...
   + y(2:n-3).^3/dx2 - c2*y(1:n-4);

% anti-diffusion terms
%     - [y(i+1) - 2*y(i) + y(i-1)]/dx2
dy(3:n-2) = -(y(4:n-1) - 2*y(3:n-2) + y(2:n-3))/dx2;


% enforce boundary conditions (u_x = u_xxx = 0)
dy(1) = 39/17*dy(3) - 28/17*dy(4) + 6/17*dy(5);  
dy(2) = 67/34*dy(3) - 21/17*dy(4) + 9/34*dy(5);  
dy(n-1) = 67/34*dy(n-2) - 21/17*dy(n-3) + 9/34*dy(n-4);  
dy(n) = 39/17*dy(n-2) - 28/17*dy(n-3) + 6/17*dy(n-4);


% end function
