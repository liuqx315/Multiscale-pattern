% driver for racing part of projective integration project
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

clear
%maxNumCompThreads(1);

% set the total integration time
Tf = 10;

% set desired output times
tout = [0,Tf];

% set the time step size bounds, tolerance
hmin = 1e-6;
hmax = 1.0;
tol = 1e-6;

% get the DIRK Butcher table
% mname = 'ARK3(2)4L[2]SA-ESDIRK';
mname = 'ARK4(3)6L[2]SA-ESDIRK';
B = butcher(mname);


printf('\nRunning tests with integrator: ')
disp(mname)
disp('   ')

%%%%%% First set of tests -- accuracy %%%%%%%

%%%%%%%%%%%%%%%%
% Test 1: initially very fast, then slower evolution
Pdata.a = 1.2; 
Pdata.b = 2.5; 
Pdata.ep = 1e-5;
u0 = 3.9;
v0 = 1.1;
w0 = 2.8;
Y0 = [u0; v0; w0];

% integrate using student's solver, and compare with "true" solution
Ytrue = load('Ytrue_test1');
err_max = 0;
err_rms = 0;
[t,Y,ns] = solve_DIRK2('fbruss1', 'Jbruss1', Pdata, tout, Y0, B, tol, hmin, hmax);
for j=1:3
   diff = (Y(j,end) - Ytrue(j))/Ytrue(j);
   err_max = max([err_max, abs(diff)]);
   err_rms = err_rms + diff^2;
end
err_rms = sqrt(err_rms/3);
fprintf('Accuracy Test 1 Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);


%%%%%%%%%%%%%%%%
% Test 2: initially slow then fast evolution
Pdata.a = 1; 
Pdata.b = 3.5; 
Pdata.ep = 5e-6;
u0 = 1.2;
v0 = 3.1;
w0 = 3;
Y0 = [u0; v0; w0];

% integrate using student's solver, and compare with "true" solution
Ytrue = load('Ytrue_test2');
err_max = 0;
err_rms = 0;
[t,Y,ns] = solve_DIRK2('fbruss2', 'Jbruss2', tout, Y0, B, tol, hmin, hmax);
for j=1:3
   diff = (Y(j,end) - Ytrue(j))/Ytrue(j);
   err_max = max([err_max, abs(diff)]);
   err_rms = err_rms + diff^2;
end
err_rms = sqrt(err_rms/3);
fprintf('Accuracy Test 2 Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);


%%%%%%%%%%%%%%%%
% Test 3: initially very rapid, then equilibrium or steady growth
Pdata.a = 0.5; 
Pdata.b = 3; 
Pdata.ep = 5e-4;
u0 = 3;
v0 = 3;
w0 = 3.5;
Y0 = [u0; v0; w0];

% integrate using student's solver, and compare with "true" solution
Ytrue = load('Ytrue_test3');
err_max = 0;
err_rms = 0;
[t,Y,ns] = solve_DIRK2('fbruss3', 'Jbruss3', tout, Y0, B, tol, hmin, hmax);
for j=1:3
   diff = (Y(j,end) - Ytrue(j))/Ytrue(j);
   err_max = max([err_max, abs(diff)]);
   err_rms = err_rms + diff^2;
end
err_rms = sqrt(err_rms/3);
fprintf('Accuracy Test 3 Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);
fprintf('\n');



% end of script
