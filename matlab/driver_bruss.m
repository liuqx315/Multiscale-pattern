% driver for racing part of projective integration project
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

clear

% set the total integration time
Tf = 10;

% set desired output times
tout = [0,Tf];

% set the time step size bounds, tolerance
hmin = 1e-6;
hmax = 1.0;
rtol = 1e-6;
atol = 1e-14*ones(3,1);

% get the DIRK Butcher tables
%mname = 'ARK3(2)4L[2]SA-ESDIRK';
%mname = 'ARK4(3)6L[2]SA-ESDIRK';
mname = 'ARK5(4)8L[2]SA-ESDIRK';
B = butcher(mname);

%mname2 = 'ARK3(2)4L[2]SA-ERK';
%mname2 = 'ARK4(3)6L[2]SA-ERK';
mname2 = 'ARK5(4)8L[2]SA-ERK';
B2 = butcher(mname2);


%%%%%% Test 1: initially very fast, then slower evolution  %%%%%%%
fprintf('\nSingle-cell brusselator test 1 (rtol = %g, atol = %g)\n',rtol,atol(1))

global Pdata;
Pdata.a = 1.2; 
Pdata.b = 2.5; 
Pdata.ep = 1e-5;
u0 = 3.9;
v0 = 1.1;
w0 = 2.8;
Y0 = [u0; v0; w0];

% get "true" solution
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s('f_bruss', tout, Y0, opts);


%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning ode15s to check accuracy/efficiency\n')

% see number of steps required by ode15s
opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
[t,Y] = ode15s('f_bruss', [0,Tf], Y0, opts);

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
[t,Y] = ode45('f_bruss', [0,Tf], Y0, opts);

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



%%%%%%%%%%%%%%%%
fprintf('\nRunning with SDIRK integrator: %s\n',mname)

% integrate using adaptive solver
[t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);

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
fprintf('\nRunning tests with ARK pair: %s / %s\n',mname,mname2)

% integrate using ARK solver
[t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
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


%%%%%%%%%%%%%%%%
fprintf('\nRunning tests with ERK integrator: %s\n',mname2)

% integrate using ERK solver
[t,Y,ns] = solve_ERK('f_bruss', 'EStab_bruss', tout, Y0, B2, rtol, atol, hmin, hmax);

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
fprintf('   work = %i\n\n',ns);



fprintf('\n---------------------------\n\n');



%%%%%% Test 2: initially slow then fast evolution  %%%%%%%
fprintf('Single-cell brusselator test 2 (rtol = %g, atol = %g)\n',rtol,atol(1))

Pdata.a = 1; 
Pdata.b = 3.5; 
Pdata.ep = 5e-6;
u0 = 1.2;
v0 = 3.1;
w0 = 3;
Y0 = [u0; v0; w0];

% get "true" solution
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s('f_bruss', tout, Y0, opts);


%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning ode15s to check accuracy/efficiency\n')

% see number of steps required by ode15s
opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
[t,Y] = ode15s('f_bruss', [0,Tf], Y0, opts);

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
[t,Y] = ode45('f_bruss', [0,Tf], Y0, opts);

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



%%%%%%%%%%%%%%%%
fprintf('\nRunning with SDIRK integrator: %s\n',mname)

% integrate using adaptive solver
[t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);

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
fprintf('\nRunning tests with ARK pair: %s / %s\n',mname,mname2)

% integrate using ARK solver
[t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', tout, Y0, B, B2, rtol, atol, hmin, hmax);

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
[t,Y,ns] = solve_ERK('f_bruss', 'EStab_bruss', tout, Y0, B2, rtol, atol, hmin, hmax);

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
fprintf('   work = %i\n\n',ns);



fprintf('\n---------------------------\n\n');



%%%%%% Test 3: initially very rapid, then equilibrium or steady growth  %%%%%%%
fprintf('Single-cell brusselator test 3 (rtol = %g, atol = %g)\n',rtol,atol(1))

Pdata.a = 0.5; 
Pdata.b = 3; 
Pdata.ep = 5e-4;
u0 = 3;
v0 = 3;
w0 = 3.5;
Y0 = [u0; v0; w0];

% get "true" solution
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s('f_bruss', tout, Y0, opts);


%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning ode15s to check accuracy/efficiency\n')

% see number of steps required by ode15s
opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
[t,Y] = ode15s('f_bruss', [0,Tf], Y0, opts);

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
[t,Y] = ode45('f_bruss', [0,Tf], Y0, opts);

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



%%%%%%%%%%%%%%%%
fprintf('\nRunning with SDIRK integrator: %s\n',mname)

% integrate using adaptive solver
[t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);

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
fprintf('\nRunning tests with ARK pair: %s / %s\n',mname,mname2)

% integrate using ARK solver
[t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
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
fprintf('   work = %i\n\n',ns);



%%%%%%%%%%%%%%%%
fprintf('\nRunning tests with ERK integrator: %s\n',mname2)

% integrate using ERK solver
[t,Y,ns] = solve_ERK('f_bruss', 'EStab_bruss', tout, Y0, B2, rtol, atol, hmin, hmax);

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
fprintf('   work = %i\n\n',ns);



% end of script
