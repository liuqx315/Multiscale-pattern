function [tvals,Y,nsteps] = solve_ARK(fcnI, fcnE, Jfcn, EStabFn, tvals, Y0, Bi, Be, rtol, atol, hmin, hmax, hmethod)
% usage: [tvals,Y,nsteps] = solve_ARK(fcnI, fcnE, Jfcn, EStabFn, tvals, Y0, Bi, Be, rtol, atol, hmin, hmax, hmethod)
%
% Adaptive time step additive RK solver for the vector-valued ODE problem
%     y' = F(t,Y), t in tspan,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:  fcnI = function name for implicit portion of ODE right-hand side
%          fcnE = function name for explicit portion of ODE right-hand side,
%                 where   F(t,Y) = fcnI(t,Y) + fcnE(t,Y)
%          Jfcn = function name for Jacobian of fcnI, J(t,Y)
%          EStabFn = function name for stability constraint on fcnE
%          tvals = [t0, t1, t2, ..., tN]
%          Y0 = initial values
%          Bi = Butcher matrix for SDIRK coefficients, of the form
%                Bi = [c A;
%                      0 b;
%                      0 b2 ]
%          Be = Butcher matrix for SDIRK coefficients, of the form
%                Be = [c A;
%                      0 b;
%                      0 b2 ]
%              The b2 row is optional in both Bi and Be, but if not provided
%              the method will default to taking fixed step sizes, at
%              min(hmin, hstab), where hstab is the result of EStabFn(t,Y).
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
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% extract DIRK method information from Bi
[Brows, Bcols] = size(Bi);
si = Bcols - 1;
ci = Bi(1:si,1);
bi = (Bi(si+1,2:si+1))';
Ai = Bi(1:si,2:si+1);
if (Brows > Bcols) 
   embedded = 1;
   bi2 = (Bi(si+2,2:si+1))';
else
   embedded = 0;
end

% extract ERK method information from Be
[Brows, Bcols2] = size(Be);
se = Bcols2 - 1;
ce = Be(1:se,1);
be = (Be(se+1,2:se+1))';
Ae = Be(1:se,2:se+1);
if ((Brows > Bcols2) && (embedded == 1))
   embedded = 1;
   be2 = (Be(se+2,2:se+1))';
else
   embedded = 0;
end

% ensure that DIRK and ERK method match at internal stage times
if (si == se-1)
   % if DIRK method not written with initial explicit stage, pad DIRK with a
   % zero first stage to match ERK
   Ai = [0, 0*bi'; 0*ci, Ai];
   ci = [0; ci];
   bi = [0, bi'];
   si = si + 1;
   if (embedded)
      bi2 = [0, bi2];
      Bi = [ci, Ai; 0, bi; 0, bi2];
   else
      Bi = [ci, Ai; 0, bi];
   end
end

% determine method order (and embedding order)
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


% ensure that ci now equals ce, otherwise return with error
if (si ~= se)
   error('incompatible SDIRK/EKR pair for ARK method, internal stage mismatch!')
end
if (norm(ci-ce) > sqrt(eps))
   error('incompatible SDIRK/EKR pair for ARK method, internal stage mismatch!')
end


% initialize output arrays
N = length(tvals);
m = length(Y0);
Y = zeros(m,N);
Y(:,1) = Y0;

% initialize diagnostics
h_s = 0;
h_a = 0;
c_fails = 0;
a_fails = 0;

% set the solver parameters
newt_maxit = 20;
%newt_tol   = 1e-10;
newt_tol   = 1e-12;
newt_alpha = 1;
dt_safety  = 0.1;
dt_reduce  = 0.1;
dt_stable  = 0.5;


% store temporary states
t = tvals(1);
Ynew = Y0;

% create Fdata structure (since used in nonlinear solver, set IRK values)
Fdata.fname = fcnI;
Fdata.Jname = Jfcn;
Fdata.B = Bi;
Fdata.s = si;

% add in explicit method data for eventual solution
Fdata.fnameE = fcnE;
Fdata.Be = Be;

% set initial time step size
h = hmin;

% initialize error weight vector
ewt = rtol*Ynew + atol;

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
      z = zeros(m,si);

      % reset step failure flag
      st_fail = 0;
      
      % loop over stages
      for stage=1:si
	 
	 % set stage initial guess as previous stage solution
	 Yguess = Ynew;
      
	 % set stage number into Fdata structure
	 Fdata.stage = stage;
	 
	 % construct RHS comprised of old time data
	 %    zi = y_n + h*sum_{j=1}^i (Ai(i,j)*fi(zj)) 
	 %             + h*sum_{j=1}^{i-1} (Ae(i,j)*fe(zj))
	 %    zi = y_n + h*Ai(i,i)*fi(zi)
	 %             + h*sum_{j=1}^{i-1} (Ai(i,j)*fi(zj) + Ae(i,j)*fe(zj))
	 %    zi - h*Ai(i,i)*fi(zi) = 
	 %             y_n + h*sum_{j=1}^{i-1} (Ai(i,j)*fi(zj) + Ae(i,j)*fe(zj))
	 % =>
	 %    rhs = y_n + h*sum_{j=1}^{i-1} (Ai(i,j)*fi(zj) + Ae(i,j)*fe(zj))
	 Fdata.rhs = Y0;
	 for j=1:stage-1
	    Fdata.rhs = Fdata.rhs + h*Ai(stage,j)*feval(fcnI,t+h*ci(j),z(:,j)) ...
		                  + h*Ae(stage,j)*feval(fcnE,t+h*ce(j),z(:,j));
	 end
	 
	 % call newton solver to compute new stage solution
	 Fdata.t = t;
	 [Ynew,inewt,ierr] = newton_damped('F_DIRK', 'A_DIRK', Yguess, ...
	     Fdata, newt_tol, newt_maxit, newt_alpha);
	 
	 % check newton error flag, if failure break out of stage loop
	 if (ierr ~= 0) 
	    % fprintf('solve_ARK warning: stage failure, cutting timestep\n');
	    st_fail = 1;
	    c_fails = c_fails + 1;
	    break;
	 end
	 
	 % update stored solution with new value
	 z(:,stage) = Ynew;
	 
      end
      nsteps = nsteps + 1;


      % compute new solution, embedding (if available)
      [Ynew,Yerr] = Y_ARK(z,Fdata);

      
      % check whether step accuracy meets tolerances (only if stages successful)
      if ((st_fail == 0) & embedded)

	 % compute error in current step
	 err_step = max(norm(Yerr./ewt,inf), eps);
	 
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
         ewt = rtol*Ynew + atol;
         
	 % for embedded methods, estimate error and update time step
	 if (embedded) 
	    h = h_estimate(Yerr, h, ewt, q_method, hmethod, 0);
	 else
	    h = hmin;
	 end
	 
	 % limit time step by explicit method stability condition
	 hstab = dt_stable * feval(EStabFn, t, Ynew);
	 if (h < hstab)
	    h_a = h_a + 1;
	 else
	    h_s = h_s + 1;
	    % fprintf('   h limited for explicit stability: h = %g, h_stab = %g\n',h,hstab);
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
fprintf('solve_ARK: took %i accuracy and %i stability steps\n',h_a,h_s);
fprintf('  step failures:  %i convergence,  %i accuracy\n',c_fails,a_fails);


% end function
