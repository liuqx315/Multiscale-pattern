function [tvals,Y,nsteps] = solve_IRK(fcn, Jfcn, tvals, Y0, B, hmax)
% usage: [tvals,Y,nsteps] = solve_IRK(fcn, Jfcn, tvals, Y0, B, hmax)
%
% IRK solver for the vector-valued ODE problem
%     y' = F(t,Y), t in tspan,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:  fcn = function name for ODE right-hand side, F(t,Y)
%          Jfcn = function name for Jacobian of ODE right-hand side, J(t,Y)
%          tvals = [t0, t1, t2, ..., tN]
%          Y0 = initial values
%          B = Butcher matrix for IRK coefficients, of the form
%                 B = [c A;
%                      0 b;
%                      0 b2 ]
%              The b2 row is optional, and provides coefficients for an
%              embedded method.
%          hmax = maximum internal time step size (can be smaller than t(i)-t(i-1))
%          nsteps = number of internal time steps taken (all stages)
%
% Outputs: t = tspan
%          y = [y(t0), y(t1), y(t2), ..., y(tN)], where each
%              y(t*) is a column vector of length m.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% get number of stages for IRK method
[Brows, Bcols] = size(B);
s = Bcols - 1;
if (Brows > Bcols) 
   embedded = 1;
else
   embedded = 0;
end

% initialize output arrays
N = length(tvals);
m = length(Y0);
Y = zeros(m,N);
Y(:,1) = Y0;

% set the solver parameters
newt_maxit = 20;
newt_tol   = 1e-10;
newt_alpha = 1;

% store temporary states
t = tvals(1);
Ynew = Y0;

% create Fdata structure
Fdata.fname = fcn;
Fdata.Jname = Jfcn;
Fdata.B = B;
Fdata.s = s;

% initialize work counter
nsteps = 0;

% iterate over time steps
for tstep = 2:length(tvals)

   % loop over internal time steps to get to desired output time
   while (t < tvals(tstep))
      
      % set internal time step
      h = min([hmax, tvals(tstep)-t]);
      Fdata.h = h;
      
      % call newton solver to update solution in time
      Fdata.yold = Ynew;
      Fdata.t = t;
      fold = feval(fcn, t, Ynew);
      z = zeros(s*m,1);
      for i = 0:s-1
	 z(i*m+1:(i+1)*m) = Ynew;
      end
      [z,ierr] = newton_damped('F_IRK', 'A_IRK', z, Fdata, newt_tol, ...
	  newt_maxit, newt_alpha);
      Ynew = Y_IRK(z,Fdata);
      
      % check newton error flag, if failure return with error
      if (ierr ~= 0) 
	 error('solve_IRK error: newton failure, reduce step and retry');
      end
	 
      % update time, work counter
      t = t + h;
      nsteps = nsteps + s;
   
   end

   % store updated solution
   Y(:,tstep) = Ynew;
   
end


% end function
