%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

clear

% set the total integration time
Tf = 10;

% set desired output times
tvals = [0,Tf];

% set the time step size bounds, tolerance
hmin = 1e-6;
hmax = 1.0;
rtol = 1e-6;
atol = 1e-14*ones(3,1);

% get the DIRK Butcher tables
%mname = 'ARK3(2)4L[2]SA-ESDIRK';
mname = 'ARK4(3)6L[2]SA-ESDIRK';
%mname = 'ARK5(4)8L[2]SA-ESDIRK';
Bi = butcher(mname);

%mname2 = 'ARK3(2)4L[2]SA-ERK';
mname2 = 'ARK4(3)6L[2]SA-ERK';
%mname2 = 'ARK5(4)8L[2]SA-ERK';
Be = butcher(mname2);


fprintf('\nSingle-cell brusselator test 1 (rtol = %g, atol = %g)\n',rtol,atol(1))
global Pdata;
Pdata.a = 1.2; 
Pdata.b = 2.5; 
Pdata.ep = 1e-5;
fcn  = 'f_bruss';
fcnI = 'fi_bruss';
fcnE = 'fe_bruss';
Jfcn = 'Ji_bruss';
EStabFn = 'EStab_bruss';
u0 = 3.9;
v0 = 1.1;
w0 = 2.8;
Y0 = [u0; v0; w0];
Yt = Y0;


%%%%%%%%%%%%%%%%
fprintf('\nRunning tests with ARK pair: %s / %s\n',mname,mname2)


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
if (se == 4) 
   q_method = 3;
   p_method = 2;
elseif (se == 6)
   q_method = 4;
   p_method = 3;
elseif (se == 8)
   q_method = 5;
   p_method = 4;
end   


% initialize output arrays
N = length(tvals);
m = length(Y0);
Y = zeros(m,N);
Y(:,1) = Y0;

% initialize diagnostics
h_e = 0;
h_i = 0;

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

      % reset stage success flag
      st_success = 0;
      
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
	    fprintf('solve_ARK warning: stage failure, cutting timestep\n');
	    st_success = 1;
	    break;
	 end
	 
	 % update stored solution with new value
	 z(:,stage) = Ynew;
	 
      end
      nsteps = nsteps + 1;

      % check whether all stages solved successfully
      if (st_success == 0) 
      
	 % use stage solutions to determine time-evolved answer, update old
	 [Ynew,Y2] = Y_ARK(z,Fdata);
	 Y0 = Ynew;

         % compute true solution at new time
         opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
         [tout,Ytrue] = ode15s(fcn, [t, t+h/2, t+h], Yt, opts);
         Yt = (Ytrue(end,:))';


	 % update time
	 t = t + h;

         % compute true and estimated errors
         e_true = max(abs((Ynew - Yt)./Yt));
         e_aprx = max(abs((Ynew - Y2)./Yt));
         fprintf(' t = %8.6f:  e_true = %7.1e,  e_aprx = %7.1e,  |fi| = %7.1e,  |fe| = %7.1e,  |f| = %7.1e\n',...
                 t, e_true, e_aprx, max(abs(feval(fcnI,t,Ynew))), ...
                 max(abs(feval(fcnE,t,Ynew))), max(abs(feval(fcn,t,Ynew))));
         
         
	 % for embedded methods, estimate error and update time step,
	 % assuming that method local truncation error order equals number of stages
	 if (embedded) 
	    h = h_estimate(Ynew, Y2, h, rtol, atol, q_method);
	 else
	    h = hmin;
	 end
	 
	 % limit time step by explicit method stability condition
	 hstab = dt_stable * feval(EStabFn, t, Ynew);
	 if (h < hstab)
	    h_i = h_i + 1;
	 else
	    h_e = h_e + 1;
	    fprintf('   h limited for explicit stability: h = %g, h_stab = %g\n',h,hstab);
	 end
	 h = min([h, hstab]);
	 
      % if any stages failed, reduce step size and retry
      else
	 
	 % reset solution guess, reduce time step
	 Ynew = Y0;
	 h = h * dt_reduce;
      
      end
      
   end

   % store updated solution
   Y(:,tstep) = Ynew;
   
end

% end of script
