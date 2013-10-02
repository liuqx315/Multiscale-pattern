function [X,Y] = stab_region(A,b,Theta)
% Usage: [X,Y] = stab_region(A,b,Theta)
% 
% Inputs:
%    A is a Butcher table matrix
%    b is a Butcher table gluing coefficients
%    Theta is an array of values in the interval [0,2*pi)
%
% We consider the RK stability function
%    R(eta) = 1 + eta*b'*((I-eta*A)\e)
%
% On the boundary of the RK stability region, |R(eta)|=1, implying
% that R(eta) = exp(i*theta), for some value of theta.
%
% For each theta in the input array Theta, and for the RK method
% defined by the array b and matrix A, we find the coordinates
% (x,y) of the complex number eta that gives rise to that point on
% the stability region boundary.  
%
% Outputs:
%    X is an array of real components of the stability boundary
%    Y is an array of imaginary components of the stability boundary
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% set perturbation, tolerance, maxit
sig = sqrt(eps);
ctol = 1e-7;
dtol = 1000;
maxit = 50;

% determine total number of plotting points
N = length(Theta);

% extract the components of the Butcher table
s = length(b);
I = eye(s);
e = ones(s,1);

% initialize the output arrays
X = zeros(N,1);
Y = zeros(N,1);

% initialize guess
x = 0;
y = 0;

% loop over plotting points
for k=1:N

   % set the root-finding and Jacobian functions for this theta
   eta = @(x,y) x+i*y;
   fc = @(x,y) 1 + eta(x,y)*b'*((I-eta(x,y)*A)\e) - exp(i*Theta(k));
   fx = @(x,y) real(fc(x,y));
   fy = @(x,y) imag(fc(x,y));
   Jxx = @(x,y) (fx(x+sig,y)-fx(x,y))/sig;
   Jxy = @(x,y) (fx(x,y+sig)-fx(x,y))/sig;
   Jyx = @(x,y) (fy(x+sig,y)-fy(x,y))/sig;
   Jyy = @(x,y) (fy(x,y+sig)-fy(x,y))/sig;
   f = @(x,y) [fx(x,y); fy(x,y)];
   J = @(x,y) [Jxx(x,y), Jxy(x,y); Jyx(x,y), Jyy(x,y)];

   % perform Newton iteration to solve for point
   for iter=1:maxit
      fun = f(x,y);
      %      fprintf('point %i, iter %i: |fun| = %g\n',k,iter,norm(fun));
      if (norm(fun) < ctol)   % check for convergence
         X(k) = x;            % store point
         Y(k) = y;
         break;               % move on to next point
      end
      if (norm(fun) > dtol)   % check for divergence
         fprintf('stab_region error: Newton divergence for point %i\n', k);
         return;              % return to calling routine
      end
      C = [x;y] - J(x,y)\fun; % perform Newton update
      x = C(1);  y = C(2);    % extract current guess
   end
   
end
   

% end of function
