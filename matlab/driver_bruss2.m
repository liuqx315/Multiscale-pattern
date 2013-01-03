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
tout = linspace(0,Tf,11);

% set the time step size bounds, tolerance
hmin = 1e-6;
hmax = 1.0;
rtol = 1e-6;
atol = 1e-10*ones(3,1);
% $$$ hmin = 0.1;
% $$$ hmax = 0.1;
% $$$ rtol = 1;
% $$$ atol = ones(3,1);
hmethod = 3;

% get the DIRK Butcher tables
mname = 'SDIRK-5-4';
B = butcher(mname);

fprintf('\nSingle-cell brusselator test 1 (rtol = %g, atol = %g)\n',rtol,atol(1))

global Pdata;
Pdata.a = 1.2; 
Pdata.b = 2.5; 
Pdata.ep = 1e-5;
%Pdata.ep = 0.01;
u0 = 3.9;
v0 = 1.1;
w0 = 2.8;
Y0 = [u0; v0; w0];

% get "true" solution
opts = odeset('RelTol',1e-12, 'AbsTol',1e-14*ones(3,1),'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s('f_bruss', tout, Y0, opts);


%%%%%%%%%%%%%%%%
fprintf('\nRunning with SDIRK integrator: %s\n',mname)

% integrate using adaptive solver
[t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, ...
                      rtol, atol, hmin, hmax, hmethod);

% compute error
err_max = max(max(abs((Y'-Ytrue)./Ytrue)));
err_rms = sum(sum(((Y'-Ytrue)./Ytrue).^2));
err_rms = sqrt(err_rms/numel(Y));
fprintf('\nResults:\n')
fprintf('      t       u        v        w      u_err      v_err      w_err\n');
fprintf('   ------------------------------------------------------------------\n');
for i=1:length(Y)
   fprintf('  %6.2f  %7.3f  %7.3f  %7.3f  %9.2e  %9.2e  %9.2e\n', ...
           tout(i), Y(1,i), Y(2,i), Y(3,i), Ytrue(i,1)-Y(1,i), ...
           Ytrue(i,2)-Y(2,i), Ytrue(i,3)-Y(3,i));
end
fprintf('   ------------------------------------------------------------------\n');

fprintf('\nOverall Accuracy/Work Results:\n')
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);



% end of script
