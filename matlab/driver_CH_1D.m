% driver for 1D Cahn-Hilliard problem,
%   u_t = -\partial_xx (c^2 \partial_xx u - u(u^2 - 1)), on [0,1]
%      u_x = 0     at x=0,x=1
%      u_xxx = 0   at x=0,x=1
%      u(x,0) = u_0 + r,  r is a uniform random number in [-0.05,0.05].
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

clear

% set the problem parameters
c = 0.00078;

% set the total integration time, mesh size
%Tf = 0.1;
Tf = 0.001;
N = 2049;

% set desired output times
tout = linspace(0,Tf,50);

% set the time step size bounds, tolerance
hmin = 1e-12;
hmax = 1e-5;
rtol = 1e-6;
atol = 1e-14*ones(N,1);

% get the DIRK Butcher tables
mname = 'ARK3(2)4L[2]SA-ESDIRK';
%mname = 'ARK4(3)6L[2]SA-ESDIRK';
B = butcher(mname);

mname2 = 'ARK3(2)4L[2]SA-ERK';
%mname2 = 'ARK4(3)6L[2]SA-ERK';
B2 = butcher(mname2);


% initial conditions
%      u(x,0) = u_0 + r,  r is a uniform random number in [-0.05,0.05].
Y0 = 0.5 + 0.01*(rand(N,1)-0.5);

% boundary conditions
%      u_x = 0     at x=0,x=1
%      u_xxx = 0   at x=0,x=1
Y0(1)   = 39/17*Y0(3)   - 28/17*Y0(4)   + 6/17*Y0(5);
Y0(2)   = 67/34*Y0(3)   - 21/17*Y0(4)   + 9/34*Y0(5);
Y0(N-1) = 67/34*Y0(N-2) - 21/17*Y0(N-3) + 9/34*Y0(N-4);
Y0(N)   = 39/17*Y0(N-2) - 28/17*Y0(N-3) + 6/17*Y0(N-4);

% store the problem parameters for function evaluations
global Pdata;
Pdata.c = c;
Pdata.dx = 1/(N-1);
Pdata.n = N;


fprintf('\n1D Cahn-Hilliard test (rtol = %g, atol = %g)\n',rtol,atol(1))

% get "true" solution
fprintf('\nComputing "true" solution with ode15s\n')
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s('f_CH_1D', tout, Y0, opts);


%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning ode15s to check accuracy/efficiency\n')

% see number of steps required by ode15s
opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
[t,Y] = ode15s('f_CH_1D', [0,Tf], Y0, opts);

% compute error
err_max = 0;
err_rms = 0;
for j=1:length(Y0)
   diff = (Y(end,j) - Ytrue(end,j))/Ytrue(end,j);
   err_max = max([err_max, abs(diff)]);
   err_rms = err_rms + diff^2;
end
err_rms = sqrt(err_rms/length(Y0));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',length(t));



%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning with SDIRK integrator: %s\n',mname)

% integrate using adaptive solver
[t,Y,ns] = solve_DIRK2('f_CH_1D', 'J_CH_1D', tout, Y0, B, rtol, atol, hmin, hmax);

% compute error
err_max = 0;
err_rms = 0;
for j=1:length(Y0)
   diff = (Y(j,end) - Ytrue(end,j))/Ytrue(end,j);
   err_max = max([err_max, abs(diff)]);
   err_rms = err_rms + diff^2;
end
err_rms = sqrt(err_rms/length(Y0));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);



%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning with ARK pair: %s / %s\n',mname,mname2)

% integrate using ARK solver
[t,Y,ns] = solve_ARK('fi_CH_1D', 'fe_CH_1D', 'Ji_CH_1D', 'EStab_CH_1D', ...
    tout, Y0, B, B2, rtol, atol, hmin, hmax);

% compute error
err_max = 0;
err_rms = 0;
for j=1:length(Y0)
   diff = (Y(j,end) - Ytrue(end,j))/Ytrue(end,j);
   err_max = max([err_max, abs(diff)]);
   err_rms = err_rms + diff^2;
end
err_rms = sqrt(err_rms/length(Y0));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);


% end of script
