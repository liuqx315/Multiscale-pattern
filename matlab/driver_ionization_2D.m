% driver for 2D radiative ionization PDE test
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

clear

% set the total integration time, mesh size
Tf = 80;
N = 1.0/0.002;

% set the brusselator parameters
alpha = 0.6; 
beta = 2; 

% set desired output times
tout = linspace(0,Tf,50);

% set the time step size bounds, tolerance
hmin = 1e-6;
hmax = 1.0;
rtol = 1e-6;
atol = 1e-14*ones(2,1);

% get the DIRK Butcher tables
% mname = 'ARK3(2)4L[2]SA-ESDIRK';
mname = 'ARK4(3)6L[2]SA-ESDIRK';
B = butcher(mname);

%mname2 = 'ARK3(2)4L[2]SA-ERK';
mname2 = 'ARK4(3)6L[2]SA-ERK';
B2 = butcher(mname2);


% initial conditions
xspan = linspace(0,1,N)';
u0 = alpha + 0.1*sin(pi*xspan);
v0 = beta/alpha + 0.1*sin(pi*xspan);
Y0 = [u0; v0];

% set problem parameter data structure
global Pdata;
Pdata.N = N;


fprintf('\n2D ionization test (rtol = %g, atol = %g)\n',rtol,atol(1))


% get "true" solution
fprintf('\nComputing "true" solution with ode15s\n')
opts = odeset('RelTol',1e-12, 'AbsTol',atol,...
              'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s('f_ionization_2D', tout, Y0, opts);


%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning ode15s to check accuracy/efficiency\n')

% see number of steps required by ode15s
opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
[t,Y] = ode15s('f_ionization_2D', [0,Tf], Y0, opts);

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
[t,Y,ns] = solve_DIRK('f_ionization_2D', 'J_ionization_2D', tout, Y0, B, ...
    rtol, atol, hmin, hmax);

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
[t,Y,ns] = solve_ARK('fi_ionization_2D', 'fe_ionization_2D', ...
    'Ji_ionization_2D', 'EStab_ionization_2D', tout, Y0, B, B2, rtol, atol, hmin, hmax);

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
