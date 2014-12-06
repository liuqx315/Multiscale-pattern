function [tvals,Y,nsteps] = solve_ERK(fcn, EStabFn, tvals, Y0, B, rtol, atol, hmin, hmax, hmethod)
% usage: [tvals,Y,nsteps] = solve_ERK(fcn, EStabFn, tvals, Y0, B, rtol, atol, hmin, hmax, hmethod)
%
% Adaptive time step explicit RK solver for the vector-valued ODE problem
%     y' = F(t,Y), t in tspan,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:  fcn = function name for ODE right-hand side,
%          EStabFn = function name for stability constraint on fcn
%          tvals = [t0, t1, t2, ..., tN]
%          Y0 = initial values
%          B = Butcher matrix for ERK coefficients, of the form
%                 B = [c A;
%                      0 b;
%                      0 b2 ]
%              The b2 row is optional, but if not provided the method will
%              default to taking fixed step sizes, at min(hmin, hstab),
%              where hstab is the result of EStabFn(t,Y). 
%          rtol = desired time accuracy relative tolerance 
%          atol = desired time accuracy absolute tolerance 
%          hmin = min internal time step size (must be smaller than t(i)-t(i-1))
%          hmax = max internal time step size (can be smaller than t(i)-t(i-1))
%          nsteps = number of internal time steps taken
%          hmethod = integer flag denoting which time adaptivity strategy to use
%
% Outputs: t = tspan
%          y = [y(t0), y(t1), y(t2), ..., y(tN)], where each
%              y(t*) is a column vector of length m.
%
% Note: to run in fixed-step mode, just call the solver using hmin=hmax as
% the desired time step size.
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% extract ERK method information from B
[Brows, Bcols] = size(B);
s = Bcols - 1;
c = B(1:s,1);
b = (B(s+1,2:s+1))';
A = B(1:s,2:s+1);
if (Brows > Bcols)
   embedded = 1;
   b2 = (B(s+2,2:s+1))';
else
   embedded = 0;
end

% estimate method order (and embedding order)
if (s == 4) 
   q_method = 3;
   p_method = 2;
elseif (s == 5)
   q_method = 4;
   p_method = 3;
elseif (s == 6)
   q_method = 4;
   p_method = 3;
elseif (s == 7)
   q_method = 5;
   p_method = 4;
elseif (s == 8)
   q_method = 5;
   p_method = 4;
else
   q_method = 5;
   p_method = 4;
end   

% initialize output arrays
N = length(tvals);
m = length(Y0);
Y = zeros(m,N);
Y(:,1) = Y0;

% initialize diagnostics
h_a = 0;
h_s = 0;
a_fails = 0;

% set the solver parameters
dt_safety  = 0.1;
dt_reduce  = 0.1;
dt_stable  = 0.5;


% store temporary states
t = tvals(1);
Ynew = Y0;

% create Fdata structure (since used in nonlinear solver, set IRK values)
Fdata.fname = fcn;
Fdata.B = B;
Fdata.s = s;

% set initial time step size
h = hmin;

% initialize error weight vector
ewt = 1.0./(rtol*Ynew + atol);

% reset time step controller
h_estimate(0, 0, 0, 0, hmethod, 1);

% initialize work counter
nsteps = 0;

% iterate over time steps
for tstep = 2:length(tvals)

   % loop over internal time steps to get to desired output time
   while (t < tvals(tstep))
      
      % bound internal time step using inputs
      h = max([h, hmin]);            % enforce minimum time step size
      h = min([h, hmax]);            % maximum time step size
      h = min([h, tvals(tstep)-t]);  % stop at output time
      Fdata.h = h;
      Fdata.yold = Y0;


      % initialize data storage for multiple stages
      z = zeros(m,s);

      % reset step failure flag
      st_fail = 0;
      
      % loop over stages
      for stage=1:s
	 
	 % set stage initial guess as previous stage solution
	 Yguess = Ynew;
      
	 % set stage number into Fdata structure
	 Fdata.stage = stage;
	 
	 % construct stage solution
	 %    zi = y_n + h*sum_{j=1}^{i-1} (A(i,j)*f(zj))
	 z(:,stage) = Y0;
	 for j=1:stage-1
	    z(:,stage) = z(:,stage) + h*A(stage,j)*feval(fcn,t+h*c(j),z(:,j));
	 end
	 
      end
      nsteps = nsteps + 1;

      % compute new solution, embedding (if available)
      Fdata.t = t;
      [Ynew,Y2] = Y_ERK(z,Fdata);

      
      % check whether step accuracy meets tolerances
      if (embedded)

	 % compute error in current step
	 err_step = max(norm((Ynew - Y2)./ewt,inf), eps);
	 
	 % if error too high, flag step as a failure (to be recomputed)
	 if (err_step > 1.2) 
	    a_fails = a_fails + 1;
	    st_fail = 1;
	 end
	 
      end
      
      % if step was successful
      if (st_fail == 0) 
      
	 % update old solution, current time
	 Y0 = Ynew;
	 t = t + h;
      
         % update error weight vector
         ewt = 1.0./(rtol*Ynew + atol);
         
	 % for embedded methods, estimate error and update time step,
	 % assuming that method local truncation error order equals number of stages
	 if (embedded) 
	    h = h_estimate(Ynew-Y2, h, ewt, q_method, hmethod, 0);
	 else
	    h = hmin;
	 end
	 
	 % limit time step by stability condition
	 hstab = dt_stable * feval(EStabFn, t, Ynew);
	 if (h < hstab)
	    h_a = h_a + 1;
	 else
	    h_s = h_s + 1;
	    %fprintf('   h limited for explicit stability: h_i = %g, h_e = %g\n',h,hstab);
	 end
	 h = min([h, hstab]);
	 
      % if step failed, reduce step size and retry
      else
	 
	 % reset solution guess, update work counter, reduce time step
	 Ynew = Y0;
	 h = h * dt_reduce;
	 h_s = h_s + 1;
         if (h <= hmin) 
            return
         end
	 
      end
      
   end

   % store updated solution
   Y(:,tstep) = Ynew;
   
end


% output diagnostics
% $$$ fprintf('solve_ERK: took %i accuracy and %i stability steps (%i failures)\n',...
% $$$     h_a, h_s, a_fails);


% end function
