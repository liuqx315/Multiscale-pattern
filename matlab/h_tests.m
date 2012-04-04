% driver that runs full suite of test problems to examine different time
% step selection strategies
clear
maxNumCompThreads(1);

% tests to run
b1_test = 1;  % brusselator test 1
b2_test = 1;  % brusselator test 2
b3_test = 1;  % brusselator test 3
v1_test = 1;  % van der Pol test 1
v2_test = 1;  % van der Pol test 2
v3_test = 1;  % van der Pol test 3

% methods to use
ode15s_meth = 1;
ode45_meth  = 0;
DIRK32_meth = 1;
DIRK43_meth = 1;
DIRK54_meth = 1;
ARK32_meth  = 1;
ARK43_meth  = 1;
ARK54_meth  = 1;
ERK32_meth  = 0;
ERK43_meth  = 0;
ERK54_meth  = 0;

% controllers to use
I_cont   = 1;   % I controller
PI_cont  = 1;   % PI controller
PID_cont = 1;   % PID controller
Gi_cont  = 1;   % Gustafsson controller
Ge_cont  = 1;   % explicit Gustafsson controller
Gie_cont = 1;   % imex Gustafsson controller

% relative tolerances to use
rtols = [1e-2, 1e-4, 1e-6];

%    start diary
!\rm -f h_test_results.txt
diary h_test_results.txt


% loop over rtols
for i_rtol = 1:length(rtols)

   % set relative tolerance for this sweep of tests
   rtol = rtols(i_rtol);

   %-----------------------------------------------------
   % brusselator test 1 -- initially fast, then slow, evolution
   if (b1_test) 
      
      % set the total integration time, desired output times, step bounds, tolerances
      Tf = 10;
      tout = linspace(0,Tf,100);
      hmin = 1e-6;
      hmax = 1.0;
      atol = 1e-14*ones(3,1);

      % set the problem parameters
      global Pdata;
      Pdata.a = 1.2; 
      Pdata.b = 2.5; 
      Pdata.ep = 1e-5;
      u0 = 3.9;
      v0 = 1.1;
      w0 = 2.8;
      Y0 = [u0; v0; w0];
      fprintf('\nSingle-cell brusselator test 1 (rtol = %g, atol = %g)\n',rtol,atol(1))

      % get "true" solution
      opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
      [t,Ytrue] = ode15s('f_bruss', tout, Y0, opts);
      if (i_rtol == 1)
	 figure()
	 plot(tout,Ytrue)
	 xlabel('t','FontSize',12), ylabel('y','FontSize',12)
	 title('Brusselator test 1','FontSize',14)
	 set(gca,'FontSize',12)
	 print('-depsc','brusselator1')
      end

      % ode15s
      if (ode15s_meth)
	 fprintf('ode15s:\n')
	 opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
	 [t,Y] = ode15s('f_bruss', [0,Tf], Y0, opts);  ns = length(t);
	 [t,Y] = ode15s('f_bruss', tout, Y0, opts);
	 fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
	     max(max(abs((Y-Ytrue)./Ytrue))),...
	     sqrt(sum(sum(((Y-Ytrue)./Ytrue).^2))/numel(Y)));
      end

      % ode45
      if (ode45_meth)
	 fprintf('ode45:\n')
	 opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
	 [t,Y] = ode45('f_bruss', [0,Tf], Y0, opts);  ns = length(t);
	 [t,Y] = ode45('f_bruss', tout, Y0, opts);
	 fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
	     max(max(abs((Y-Ytrue)./Ytrue))),...
	     sqrt(sum(sum(((Y-Ytrue)./Ytrue).^2))/numel(Y)));
      end

      if (DIRK32_meth)
	 mname = 'ARK3(2)4L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK3(2), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK3(2), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK3(2), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK32 method

      
      if (DIRK43_meth)
	 mname = 'ARK4(3)6L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK4(3), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK4(3), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK4(3), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK4(3) method

      
      if (DIRK54_meth)
	 mname = 'ARK5(4)8L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK5(4), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK5(4), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK5(4), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK5(4) method


      if (ARK32_meth)
	 mname = 'ARK3(2)4L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK3(2)4L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK3(2), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK3(2), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK3(2), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK3(2), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK3(2), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK3(2) method


      if (ARK43_meth)
	 mname = 'ARK4(3)6L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK4(3)6L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK4(3), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK4(3), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK4(3), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK4(3), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK4(3), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK4(3) method


      if (ARK54_meth)
	 mname = 'ARK5(4)8L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK5(4)8L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK5(4), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK5(4), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK5(4), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK5(4), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK5(4), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK5(4) method


      if (ERK32_meth)
	 mname2 = 'ARK3(2)4L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK3(2), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK3(2), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK3(2), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK3(2) method


      if (ERK43_meth)
	 mname2 = 'ARK4(3)6L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK4(3), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK4(3), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK4(3), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK4(3) method


      if (ERK54_meth)
	 mname2 = 'ARK5(4)8L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK5(4), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK5(4), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK5(4), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK5(4) method

      fprintf('\n\n---------------------------\n\n');

   end   % brusselator test 1

   %-----------------------------------------------------
   % brusselator test 2 -- initially slow, then fast, evolution
   if (b2_test) 

      % set the total integration time, desired output times, step bounds, tolerances
      Tf = 10;
      tout = linspace(0,Tf,100);
      hmin = 1e-6;
      hmax = 1.0;
      atol = 1e-14*ones(3,1);

      % set the problem parameters
      global Pdata;
      Pdata.a = 1; 
      Pdata.b = 3.5; 
      Pdata.ep = 5e-6;
      u0 = 1.2;
      v0 = 3.1;
      w0 = 3;
      Y0 = [u0; v0; w0];
      fprintf('\nSingle-cell brusselator test 2 (rtol = %g, atol = %g)\n',rtol,atol(1))

      % get "true" solution
      opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
      [t,Ytrue] = ode15s('f_bruss', tout, Y0, opts);
      if (i_rtol == 1)
	 figure()
	 plot(tout,Ytrue)
	 xlabel('t','FontSize',12), ylabel('y','FontSize',12)
	 title('Brusselator test 2','FontSize',14)
	 set(gca,'FontSize',12)
	 print('-depsc','brusselator2')
      end

      % ode15s
      if (ode15s_meth)
	 fprintf('ode15s:\n')
	 opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
	 [t,Y] = ode15s('f_bruss', [0,Tf], Y0, opts);  ns = length(t);
	 [t,Y] = ode15s('f_bruss', tout, Y0, opts);
	 fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
	     max(max(abs((Y-Ytrue)./Ytrue))),...
	     sqrt(sum(sum(((Y-Ytrue)./Ytrue).^2))/numel(Y)));
      end

      % ode45
      if (ode45_meth)
	 fprintf('ode45:\n')
	 opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
	 [t,Y] = ode45('f_bruss', [0,Tf], Y0, opts);  ns = length(t);
	 [t,Y] = ode45('f_bruss', tout, Y0, opts);
	 fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
	     max(max(abs((Y-Ytrue)./Ytrue))),...
	     sqrt(sum(sum(((Y-Ytrue)./Ytrue).^2))/numel(Y)));
      end


      if (DIRK32_meth)
	 mname = 'ARK3(2)4L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK3(2), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK3(2), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK3(2), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK32 method


      if (DIRK43_meth)
	 mname = 'ARK4(3)6L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK4(3), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK4(3), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK4(3), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK4(3) method


      if (DIRK54_meth)
	 mname = 'ARK5(4)8L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK5(4), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK5(4), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK5(4), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK5(4) method


      if (ARK32_meth)
	 mname = 'ARK3(2)4L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK3(2)4L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK3(2), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK3(2), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK3(2), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK3(2), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK3(2), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK3(2) method


      if (ARK43_meth)
	 mname = 'ARK4(3)6L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK4(3)6L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK4(3), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK4(3), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK4(3), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK4(3), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK4(3), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK4(3) method


      if (ARK54_meth)
	 mname = 'ARK5(4)8L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK5(4)8L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK5(4), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK5(4), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK5(4), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK5(4), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK5(4), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK5(4) method


      if (ERK32_meth)
	 mname2 = 'ARK3(2)4L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK3(2), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK3(2), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK3(2), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK3(2) method


      if (ERK43_meth)
	 mname2 = 'ARK4(3)6L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK4(3), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK4(3), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK4(3), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK4(3) method


      if (ERK54_meth)
	 mname2 = 'ARK5(4)8L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK5(4), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK5(4), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK5(4), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK5(4) method

      fprintf('\n\n---------------------------\n\n');

   end   % brusselator test 2


   %-----------------------------------------------------
   % brusselator test 3 -- initially rapid, then equilibrium or steady growth
   if (b3_test) 

      % set the total integration time, desired output times, step bounds, tolerances
      Tf = 10;
      tout = linspace(0,Tf,100);
      hmin = 1e-6;
      hmax = 1.0;
      atol = 1e-14*ones(3,1);

      % set the problem parameters
      global Pdata;
      Pdata.a = 0.5; 
      Pdata.b = 3; 
      Pdata.ep = 5e-4;
      u0 = 3;
      v0 = 3;
      w0 = 3.5;
      Y0 = [u0; v0; w0];
      fprintf('\nSingle-cell brusselator test 3 (rtol = %g, atol = %g)\n',rtol,atol(1))

      % get "true" solution
      opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
      [t,Ytrue] = ode15s('f_bruss', tout, Y0, opts);
      if (i_rtol == 1)
	 figure()
	 plot(tout,Ytrue)
	 xlabel('t','FontSize',12), ylabel('y','FontSize',12)
	 title('Brusselator test 3','FontSize',14)
	 set(gca,'FontSize',12)
	 print('-depsc','brusselator3')
      end

      % ode15s
      if (ode15s_meth)
	 fprintf('ode15s:\n')
	 opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
	 [t,Y] = ode15s('f_bruss', [0,Tf], Y0, opts);  ns = length(t);
	 [t,Y] = ode15s('f_bruss', tout, Y0, opts);
	 fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
	     max(max(abs((Y-Ytrue)./Ytrue))),...
	     sqrt(sum(sum(((Y-Ytrue)./Ytrue).^2))/numel(Y)));
      end

      % ode45
      if (ode45_meth)
	 fprintf('ode45:\n')
	 opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
	 [t,Y] = ode45('f_bruss', [0,Tf], Y0, opts);  ns = length(t);
	 [t,Y] = ode45('f_bruss', tout, Y0, opts);
	 fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
	     max(max(abs((Y-Ytrue)./Ytrue))),...
	     sqrt(sum(sum(((Y-Ytrue)./Ytrue).^2))/numel(Y)));
      end


      if (DIRK32_meth)
	 mname = 'ARK3(2)4L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK3(2), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK3(2), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK3(2), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK32 method


      if (DIRK43_meth)
	 mname = 'ARK4(3)6L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK4(3), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK4(3), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK4(3), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK4(3) method


      if (DIRK54_meth)
	 mname = 'ARK5(4)8L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK5(4), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK5(4), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK5(4), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_bruss', 'J_bruss', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK5(4) method


      if (ARK32_meth)
	 mname = 'ARK3(2)4L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK3(2)4L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK3(2), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK3(2), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK3(2), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK3(2), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK3(2), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK3(2) method


      if (ARK43_meth)
	 mname = 'ARK4(3)6L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK4(3)6L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK4(3), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK4(3), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK4(3), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK4(3), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK4(3), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK4(3) method


      if (ARK54_meth)
	 mname = 'ARK5(4)8L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK5(4)8L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK5(4), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK5(4), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK5(4), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK5(4), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK5(4), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_bruss', 'fe_bruss', 'Ji_bruss', 'EStab_bruss', ...
		tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK5(4) method


      if (ERK32_meth)
	 mname2 = 'ARK3(2)4L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK3(2), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK3(2), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK3(2), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK3(2) method


      if (ERK43_meth)
	 mname2 = 'ARK4(3)6L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK4(3), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK4(3), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK4(3), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK4(3) method


      if (ERK54_meth)
	 mname2 = 'ARK5(4)8L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK5(4), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK5(4), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK5(4), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_bruss', 'EStab_VanderPol', tout, Y0, ...
		B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)./Ytrue))),...
		sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK5(4) method

      fprintf('\n\n---------------------------\n\n');

   end   % brusselator test 3


   %-----------------------------------------------------
   % Van der Pol test 1
   if (v1_test) 

      % set the total integration time, desired output times, step bounds, tolerances
      Tf = 1.5;
      tout = linspace(0,Tf,100);
      hmin = 1e-7;
      hmax = 0.1;
      atol = 1e-14*ones(2,1);

      % set the problem parameters
      global Pdata;
      Pdata.ep = 1e-3;
      Pdata.ep2 = 1.0;
      u0 = -2;
      v0 = -2.3553013976081;
      Y0 = [u0; v0];
      fprintf('\nVan der Pol test 1 (rtol = %g, atol = %g)\n',rtol,atol(1))

      % get "true" solution
      opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
      [t,Ytrue] = ode15s('f_VanderPol', tout, Y0, opts);
      if (i_rtol == 1)
	 figure()
	 plot(tout,Ytrue)
	 xlabel('t','FontSize',12), ylabel('y','FontSize',12)
	 title('Van der Pol test 1','FontSize',14)
	 set(gca,'FontSize',12)
	 print('-depsc','vanderPol1')
      end

      % ode15s
      if (ode15s_meth)
	 fprintf('ode15s:\n')
	 opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
	 [t,Y] = ode15s('f_VanderPol', [0,Tf], Y0, opts);  ns = length(t);
	 [t,Y] = ode15s('f_VanderPol', tout, Y0, opts);
	 fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
	     max(max(abs((Y-Ytrue)))),...
	     sqrt(sum(sum(((Y-Ytrue)).^2))/numel(Y)));
      end

      % ode45
      if (ode45_meth)
	 fprintf('ode45:\n')
	 opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
	 [t,Y] = ode45('f_VanderPol', [0,Tf], Y0, opts);  ns = length(t);
	 [t,Y] = ode45('f_VanderPol', tout, Y0, opts);
	 fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
	     max(max(abs((Y-Ytrue)))),...
	     sqrt(sum(sum(((Y-Ytrue)).^2))/numel(Y)));
      end


      if (DIRK32_meth)
	 mname = 'ARK3(2)4L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK3(2), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK3(2), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK3(2), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK32 method


      if (DIRK43_meth)
	 mname = 'ARK4(3)6L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK4(3), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK4(3), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK4(3), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK4(3) method


      if (DIRK54_meth)
	 mname = 'ARK5(4)8L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK5(4), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK5(4), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK5(4), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK5(4) method


      if (ARK32_meth)
	 mname = 'ARK3(2)4L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK3(2)4L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK3(2), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK3(2), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK3(2), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK3(2), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK3(2), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK3(2) method


      if (ARK43_meth)
	 mname = 'ARK4(3)6L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK4(3)6L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK4(3), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK4(3), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK4(3), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK4(3), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK4(3), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK4(3) method


      if (ARK54_meth)
	 mname = 'ARK5(4)8L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK5(4)8L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK5(4), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK5(4), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK5(4), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK5(4), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK5(4), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK5(4) method


      if (ERK32_meth)
	 mname2 = 'ARK3(2)4L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK3(2), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK3(2), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK3(2), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK3(2) method


      if (ERK43_meth)
	 mname2 = 'ARK4(3)6L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK4(3), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK4(3), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK4(3), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK4(3) method


      if (ERK54_meth)
	 mname2 = 'ARK5(4)8L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK5(4), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK5(4), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK5(4), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK5(4) method

      fprintf('\n\n---------------------------\n\n');

   end   % van der Pol test 1


   %-----------------------------------------------------
   % Van der Pol test 2
   if (v2_test) 

      % set the total integration time, desired output times, step bounds, tolerances
      Tf = 12;
      tout = linspace(0,Tf,100);
      hmin = 1e-6;
      hmax = 1.0;
      atol = 1e-14*ones(2,1);

      % set the problem parameters
      global Pdata;
      Pdata.ep = 0.2;
      Pdata.ep2 = 0.2;
      u0 = 2;
      v0 = 0;
      Y0 = [u0; v0];
      fprintf('\nVan der Pol test 2 (rtol = %g, atol = %g)\n',rtol,atol(1))

      % get "true" solution
      opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
      [t,Ytrue] = ode15s('f_VanderPol', tout, Y0, opts);
      if (i_rtol == 1)
	 figure()
	 plot(tout,Ytrue)
	 xlabel('t','FontSize',12), ylabel('y','FontSize',12)
	 title('Van der Pol test 2','FontSize',14)
	 set(gca,'FontSize',12)
	 print('-depsc','vanderPol2')
      end

      % ode15s
      if (ode15s_meth)
	 fprintf('ode15s:\n')
	 opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
	 [t,Y] = ode15s('f_VanderPol', [0,Tf], Y0, opts);  ns = length(t);
	 [t,Y] = ode15s('f_VanderPol', tout, Y0, opts);
	 fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
	     max(max(abs((Y-Ytrue)))),...
	     sqrt(sum(sum(((Y-Ytrue)).^2))/numel(Y)));
      end

      % ode45
      if (ode45_meth)
	 fprintf('ode45:\n')
	 opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
	 [t,Y] = ode45('f_VanderPol', [0,Tf], Y0, opts);  ns = length(t);
	 [t,Y] = ode45('f_VanderPol', tout, Y0, opts);
	 fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
	     max(max(abs((Y-Ytrue)))),...
	     sqrt(sum(sum(((Y-Ytrue)).^2))/numel(Y)));
      end


      if (DIRK32_meth)
	 mname = 'ARK3(2)4L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK3(2), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK3(2), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK3(2), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK32 method


      if (DIRK43_meth)
	 mname = 'ARK4(3)6L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK4(3), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK4(3), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK4(3), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK4(3) method


      if (DIRK54_meth)
	 mname = 'ARK5(4)8L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK5(4), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK5(4), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK5(4), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK5(4) method


      if (ARK32_meth)
	 mname = 'ARK3(2)4L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK3(2)4L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK3(2), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK3(2), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK3(2), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK3(2), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK3(2), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK3(2) method


      if (ARK43_meth)
	 mname = 'ARK4(3)6L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK4(3)6L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK4(3), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK4(3), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK4(3), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK4(3), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK4(3), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK4(3) method


      if (ARK54_meth)
	 mname = 'ARK5(4)8L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK5(4)8L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK5(4), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK5(4), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK5(4), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK5(4), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK5(4), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK5(4) method


      if (ERK32_meth)
	 mname2 = 'ARK3(2)4L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK3(2), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK3(2), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK3(2), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK3(2) method


      if (ERK43_meth)
	 mname2 = 'ARK4(3)6L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK4(3), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK4(3), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK4(3), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK4(3) method


      if (ERK54_meth)
	 mname2 = 'ARK5(4)8L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK5(4), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK5(4), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK5(4), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK5(4) method

      fprintf('\n\n---------------------------\n\n');

   end   % van der Pol test 2


   %-----------------------------------------------------
   % Van der Pol test 3
   if (v3_test) 

      % set the total integration time, desired output times, step bounds, tolerances
      Tf = 1000;
      tout = linspace(0,Tf,100);
      hmin = 1e-6;
      hmax = 10;
      atol = 1e-14*ones(2,1);

      % set the problem parameters
      global Pdata;
      Pdata.ep = 1e-3;
      Pdata.ep2 = 1e-3;
      u0 = 2;
      v0 = 0;
      Y0 = [u0; v0];
      fprintf('\nVan der Pol test 3 (rtol = %g, atol = %g)\n',rtol,atol(1))

      % get "true" solution
      opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
      [t,Ytrue] = ode15s('f_VanderPol', tout, Y0, opts);
      if (i_rtol == 1)
	 figure()
	 plot(tout,Ytrue)
	 xlabel('t','FontSize',12), ylabel('y','FontSize',12)
	 title('Van der Pol test 3','FontSize',14)
	 set(gca,'FontSize',12)
	 print('-depsc','vanderPol3')
      end

      % ode15s
      if (ode15s_meth)
	 fprintf('ode15s:\n')
	 opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
	 [t,Y] = ode15s('f_VanderPol', [0,Tf], Y0, opts);  ns = length(t);
	 [t,Y] = ode15s('f_VanderPol', tout, Y0, opts);
	 fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
	     max(max(abs((Y-Ytrue)))),...
	     sqrt(sum(sum(((Y-Ytrue)).^2))/numel(Y)));
      end

      % ode45
      if (ode45_meth)
	 fprintf('ode45:\n')
	 opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
	 [t,Y] = ode45('f_VanderPol', [0,Tf], Y0, opts);  ns = length(t);
	 [t,Y] = ode45('f_VanderPol', tout, Y0, opts);
	 fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
	     max(max(abs((Y-Ytrue)))),...
	     sqrt(sum(sum(((Y-Ytrue)).^2))/numel(Y)));
      end


      if (DIRK32_meth)
	 mname = 'ARK3(2)4L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK3(2), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK3(2), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK3(2), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK32 method


      if (DIRK43_meth)
	 mname = 'ARK4(3)6L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK4(3), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK4(3), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK4(3), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK4(3) method


      if (DIRK54_meth)
	 mname = 'ARK5(4)8L[2]SA-ESDIRK';
	 B = butcher(mname);

	 if (I_cont)
	    % SDIRK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('SDIRK5(4), I controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % SDIRK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('SDIRK5(4), PI controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % SDIRK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('SDIRK5(4), PID controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % SDIRK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('SDIRK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % SDIRK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('SDIRK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_DIRK('f_VanderPol', 'J_VanderPol', tout, Y0, B, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % SDIRK5(4) method


      if (ARK32_meth)
	 mname = 'ARK3(2)4L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK3(2)4L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK3(2), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK3(2), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK3(2), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK3(2), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK3(2), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK3(2) method


      if (ARK43_meth)
	 mname = 'ARK4(3)6L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK4(3)6L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK4(3), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK4(3), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK4(3), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK4(3), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK4(3), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK4(3) method


      if (ARK54_meth)
	 mname = 'ARK5(4)8L[2]SA-ESDIRK';
	 B = butcher(mname);
	 mname2 = 'ARK5(4)8L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ARK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ARK5(4), I controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ARK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ARK5(4), PI controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ARK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ARK5(4), PID controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ARK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ARK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ARK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ARK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont

	 if (Gie_cont)
	    % ARK5(4), imex Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_imex.m h_estimate.m
	    fprintf('ARK5(4), imex Gustafsson controller:\n')
	    [t,Y,ns] = solve_ARK('fi_VanderPol', 'fe_VanderPol', 'Ji_VanderPol', ...
		'EStab_VanderPol', tout, Y0, B, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gie_cont
      end  % ARK5(4) method


      if (ERK32_meth)
	 mname2 = 'ARK3(2)4L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK3(2), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK3(2), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK3(2), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK3(2), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK3(2), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK3(2), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK3(2), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK3(2), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK3(2), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK3(2), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK3(2) method


      if (ERK43_meth)
	 mname2 = 'ARK4(3)6L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK4(3), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK4(3), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK4(3), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK4(3), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK4(3), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK4(3), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK4(3), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK4(3), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK4(3), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK4(3), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK4(3) method


      if (ERK54_meth)
	 mname2 = 'ARK5(4)8L[2]SA-ERK';
	 B2 = butcher(mname2);

	 if (I_cont)
	    % ERK5(4), I controller
	    !ln -fs h_estimate_I.m h_estimate.m
	    fprintf('ERK5(4), I controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % I_cont

	 if (PI_cont)
	    % ERK5(4), PI controller
	    !ln -fs h_estimate_PI.m h_estimate.m
	    fprintf('ERK5(4), PI controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PI_cont

	 if (PID_cont)
	    % ERK5(4), PID controller
	    !ln -fs h_estimate_PID.m h_estimate.m
	    fprintf('ERK5(4), PID controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % PID_cont

	 if (Gi_cont)
	    % ERK5(4), Gustafsson controller
	    !ln -fs h_estimate_Gustafsson.m h_estimate.m
	    fprintf('ERK5(4), Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Gi_cont

	 if (Ge_cont)
	    % ERK5(4), explicit Gustafsson controller
	    !ln -fs h_estimate_Gustafsson_exp.m h_estimate.m
	    fprintf('ERK5(4), explicit Gustafsson controller:\n')
	    [t,Y,ns] = solve_ERK('f_VanderPol', 'EStab_VanderPol', tout, ...
		Y0, B2, rtol, atol, hmin, hmax);
	    fprintf('   work = %i,   maxerr = %.2e,   rmserr = %.2e\n',ns,...
		max(max(abs((Y'-Ytrue)))),...
		sqrt(sum(sum(((Y'-Ytrue)).^2))/numel(Y)));
	 end  % Ge_cont
      end  % ERK5(4) method

      fprintf('\n\n---------------------------\n\n');

   end   % van der Pol test 2


end  % i_rtol loop


% stop diary
disp('   ')
diary off


% end of testing script