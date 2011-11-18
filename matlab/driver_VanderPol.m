% driver for stiff ODE version of Van der Pol test problem:
%    u' = v
%    v' = (v - v*u^2 - u)/ep
% where u(0) = 2,  v(0) = -0.6666654321121172, and ep = 1e-5.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

clear

% set problem parameters
ep = 1e-4;

% set the total integration time
Tf = 1.5;

% set desired output times
tout = [0,Tf];

% set the time step size bounds, tolerance
%hmin = 10^(-4.5);
hmin = 1e-6;
hmax = 0.1;
rtol = 1e-6;
atol = 1e-14*ones(2,1);

% get the DIRK Butcher tables
% mname = 'ARK3(2)4L[2]SA-ESDIRK';
mname = 'ARK4(3)6L[2]SA-ESDIRK';
B = butcher(mname);

% mname2 = 'ARK3(2)4L[2]SA-ERK';
mname2 = 'ARK4(3)6L[2]SA-ERK';
B2 = butcher(mname2);


% set the initial conditions
u0 = -2;
v0 = -2.3553013976081;
Y0 = [u0; v0];

% store problem parameters
global Pdata;
Pdata.ep = ep;


fprintf('\nVan der Pol test (rtol = %g, atol = %g)\n',rtol,atol(1))

% get "true" solution
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s('f_VanderPol', tout, Y0, opts);


%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning ode15s to check accuracy/efficiency\n')

% see number of steps required by ode15s
opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
[t,Y] = ode15s('f_VanderPol', [0,Tf], Y0, opts);

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
fprintf('\nRunning ode45 to check accuracy/efficiency\n')

% see number of steps required by ode15s
opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
[t,Y] = ode45('f_VanderPol', [0,Tf], Y0, opts);

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
[t,Y,ns] = solve_DIRK2('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, ...
    atol, hmin, hmax);

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
[t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
    'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);

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


%%%%%%%%%%%%%%%%
fprintf('\nRunning tests with ERK integrator: %s\n',mname2)

% integrate using ERK solver
[t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, Y0, B2, rtol, atol, hmin, hmax);

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
