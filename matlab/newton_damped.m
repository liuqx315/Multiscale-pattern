function [y,inewt,ierr] = newton_damped(Fcn, Afn, y0, Fdata, tol, maxit, alpha)
% usage: [y,inewt,ierr] = newton_damped(Fcn, Afn, y0, Fdata, tol, maxit, alpha)
%
% Damped Newton iteration for the nonlinear function defined by the
% function Fcn,
%     F(y,Fdata) = 0
%
% Inputs:  Fcn = function name for nonlinear residual, F(y,Fdata).  Note, all 
%                data required to evaluate F (other than y) should be stored
%                in the data structure Fdata.
%          Afn = function name for Jacobian of nonlinear residual, 
%                A = partial_y F(y,Fdata).  
%                Afn should use the same data structure for additional data
%                as F.
%          y0 = initial guess
%          Fdata = structure containing extra information for evaluating F.
%          tol = desired nonlinear tolerance
%          maxit = maximum allowed iterations
%          alpha = damping parameter, 0 < alpha <= 1
% Outputs: y = solution to root-finding problem
%          inewt = number of Newton iterations to convergence
%          ierr = output flag denoting success (0) or failure (1)
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

DO_OUTPUT = 0;

% check solver inputs
if (maxit < 1) 
   error('newton_damped error: requires at least 1 iteration (maxit)');
end
if (tol <= 0) 
   error('newton_damped error: tolerance must be positive (tol)');
end
if ((alpha <= 0) || (alpha > 1)) 
   error('newton_damped error: damping parameter not in (0,1] (alpha)');
end

% initialize result, increment vector
y = y0;
s = ones(size(y));

% perform iterations
for inewt=1:maxit
   
   % check residual for stopping
   if (norm(s,inf) < tol)
%      fprintf('  newton_damped: converged to tol %g in %i iters\n',norm(s,inf),i-1);
      ierr = 0;
      return
   end

   % compute residual, Jacobian
   F = feval(Fcn,y,Fdata);
   A = feval(Afn,y,Fdata);
   
   if (DO_OUTPUT)
      fprintf('       newt rhs = ');
      for entry=1:length(F), fprintf('%19.16g, ',F(entry)); end
      fprintf('\n');
   end

   % perform Newton update
   s = A\F;
   y = y - alpha*s;
   
end

% if we've made it to this point, the Newton iteration did not converge
ierr = 1;
%fprintf('\nnewton_damped warning: nonconvergence after %i iterations (|F| = %g)\n',maxit,norm(F,inf));

% end of function
