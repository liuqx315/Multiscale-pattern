%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------
% driver for stiff ODE version of analytical test problem:
%    u' = lam*u + 1/(1+t^2) - lam*atan(t)
% where u(0) = 0, lam = -1/ep, and ep = 1e-5.

clear

% set problem parameters
ep = 1e-5;

% set the total integration time
Tf = 10;

% set desired output times
tout = linspace(0,Tf,100);

% set the time step size bounds, tolerance
hmin = 1e-6;
hmax = 1.0;
rtol = 1e-6;
atol = 1e-14*ones(1,1);

% get the DIRK Butcher tables
mname = 'ARK3(2)4L[2]SA-ESDIRK';
%mname = 'ARK4(3)6L[2]SA-ESDIRK';
%mname = 'ARK5(4)8L[2]SA-ESDIRK';
B = butcher(mname);

mname2 = 'ARK3(2)4L[2]SA-ERK';
%mname2 = 'ARK4(3)6L[2]SA-ERK';
%mname2 = 'ARK5(4)8L[2]SA-ERK';
B2 = butcher(mname2);


% set the initial conditions
Y0 = 0;

% store problem parameters
global Pdata;
Pdata.ep = ep;


fprintf('\nAnalytical test (rtol = %g, atol = %g)\n',rtol,atol(1))

% get "true" solution
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s('f_analytic', tout, Y0, opts);


%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning ode15s to check accuracy/efficiency\n')

% see number of steps required by ode15s
opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
[t,Y] = ode15s('f_analytic', [0,Tf], Y0, opts);  nsteps = length(t);
[t,Y] = ode15s('f_analytic', tout, Y0, opts);

% compute error
err_max = max(max(abs((Y-Ytrue)./Ytrue)));
err_rms = sum(sum(((Y-Ytrue)./Ytrue).^2));
err_rms = sqrt(err_rms/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',nsteps);



%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning ode45 to check accuracy/efficiency\n')

% see number of steps required by ode15s
opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
[t,Y] = ode45('f_analytic', [0,Tf], Y0, opts);  nsteps = length(t);
[t,Y] = ode45('f_analytic', tout, Y0, opts);

% compute error
err_max = max(max(abs((Y-Ytrue)./Ytrue)));
err_rms = sum(sum(((Y-Ytrue)./Ytrue).^2));
err_rms = sqrt(err_rms/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',nsteps);



%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning with SDIRK integrator: %s\n',mname)

% integrate using adaptive solver
[t,Y,ns] = solve_DIRK('f_analytic', 'J_analytic', tout, Y0, B, rtol, ...
    atol, hmin, hmax);

% compute error
err_max = max(max(abs((Y'-Ytrue)./Ytrue)));
err_rms = sum(sum(((Y'-Ytrue)./Ytrue).^2));
err_rms = sqrt(err_rms/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);


%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning with ARK pair: %s / %s\n',mname,mname2)

% integrate using ARK solver
[t,Y,ns] = solve_ARK('fi_analytic', 'fe_analytic', 'Ji_analytic', ...
    'EStab_analytic', tout, Y0, B, B2, rtol, atol, hmin, hmax);

% compute error
err_max = max(max(abs((Y'-Ytrue)./Ytrue)));
err_rms = sum(sum(((Y'-Ytrue)./Ytrue).^2));
err_rms = sqrt(err_rms/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);


%%%%%%%%%%%%%%%%
fprintf('\nRunning tests with ERK integrator: %s\n',mname2)

% integrate using ERK solver
[t,Y,ns] = solve_ERK('f_analytic', 'EStab_analytic', tout, Y0, B2, rtol, atol, hmin, hmax);

% compute error
err_max = max(max(abs((Y'-Ytrue)./Ytrue)));
err_rms = sum(sum(((Y'-Ytrue)./Ytrue).^2));
err_rms = sqrt(err_rms/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);



% end of script
