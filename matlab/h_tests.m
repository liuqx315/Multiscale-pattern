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
tests = [b1_test, b2_test, b3_test, v1_test, v2_test, v3_test];

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

% flag to generate plots
show_plots = 0

% relative tolerances to use
rtols = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6];

%    start diary
!\rm -f h_test_results.txt
diary h_test_results.txt


% loop over test problems
for i_test = 1:length(tests)

   % if problem is active, set parameters and functions, and proceed
   if (tests(i_test) > 0) 
      
      % brusselator test 1
      switch (i_test)
         
        case(1)  % brusselator 1 -- initially fast, then slow, evolution
          fprintf('\nSingle-cell brusselator test 1\n')
          fn = 'f_bruss';
          Jn = 'J_bruss';
          fi = 'fi_bruss';
          fe = 'fe_bruss';
          Ji = 'Ji_bruss';
          Es = 'EStab_bruss';
          Tf = 10;
          tout = linspace(0,Tf,100);
          hmin = 1e-6;
          hmax = 1.0;
          atol = 1e-14*ones(3,1);
          global Pdata;
          Pdata.a = 1.2; 
          Pdata.b = 2.5; 
          Pdata.ep = 1e-5;
          u0 = 3.9;
          v0 = 1.1;
          w0 = 2.8;
          Y0 = [u0; v0; w0];

          % plot "true" solution
          if (show_plots == 1)
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
             [t,Ytrue] = ode15s('f_bruss', tout, Y0, opts);
             figure()
             plot(tout,Ytrue)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Brusselator test 1','FontSize',14)
             set(gca,'FontSize',12)
             print('-depsc','brusselator1')
          end
          

        case(2)  % brusselator 2
          fprintf('\nSingle-cell brusselator test 2\n')
          fn = 'f_bruss';
          Jn = 'J_bruss';
          fi = 'fi_bruss';
          fe = 'fe_bruss';
          Ji = 'Ji_bruss';
          Es = 'EStab_bruss';
          Tf = 10;
          tout = linspace(0,Tf,100);
          hmin = 1e-6;
          hmax = 1.0;
          atol = 1e-14*ones(3,1);
          global Pdata;
          Pdata.a = 1; 
          Pdata.b = 3.5; 
          Pdata.ep = 5e-6;
          u0 = 1.2;
          v0 = 3.1;
          w0 = 3;
          Y0 = [u0; v0; w0];

          % plot "true" solution
          if (show_plots == 1)
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
             [t,Ytrue] = ode15s('f_bruss', tout, Y0, opts);
             figure()
             plot(tout,Ytrue)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Brusselator test 2','FontSize',14)
             set(gca,'FontSize',12)
             print('-depsc','brusselator2')
          end


        case(3)  % brusselator 3
          fprintf('\nSingle-cell brusselator test 3\n')
          fn = 'f_bruss';
          Jn = 'J_bruss';
          fi = 'fi_bruss';
          fe = 'fe_bruss';
          Ji = 'Ji_bruss';
          Es = 'EStab_bruss';
          Tf = 10;
          tout = linspace(0,Tf,100);
          hmin = 1e-6;
          hmax = 1.0;
          atol = 1e-14*ones(3,1);
          global Pdata;
          Pdata.a = 0.5; 
          Pdata.b = 3; 
          Pdata.ep = 5e-4;
          u0 = 3;
          v0 = 3;
          w0 = 3.5;
          Y0 = [u0; v0; w0];

          % plot "true" solution
          if (show_plots == 1)
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
             [t,Ytrue] = ode15s('f_bruss', tout, Y0, opts);
             figure()
             plot(tout,Ytrue)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Brusselator test 3','FontSize',14)
             set(gca,'FontSize',12)
             print('-depsc','brusselator3')
          end
        
        
        case(4)  % van der Pol 1
          fprintf('\nVan der Pol test 1\n')
          fn = 'f_VanderPol';
          Jn = 'J_VanderPol';
          fi = 'fi_VanderPol';
          fe = 'fe_VanderPol';
          Ji = 'Ji_VanderPol';
          Es = 'EStab_VanderPol';
          Tf = 1.5;
          tout = linspace(0,Tf,100);
          hmin = 1e-7;
          hmax = 0.1;
          atol = 1e-14*ones(2,1);
          global Pdata;
          Pdata.ep = 1e-3;
          Pdata.ep2 = 1.0;
          u0 = -2;
          v0 = -2.3553013976081;
          Y0 = [u0; v0];

          % plot "true" solution
          if (show_plots == 1)
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
             [t,Ytrue] = ode15s('f_VanderPol', tout, Y0, opts);
             figure()
             plot(tout,Ytrue)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Van der Pol test 1','FontSize',14)
             set(gca,'FontSize',12)
             print('-depsc','vanderPol1')
          end

        
        case(5)  % van der Pol 2
          fprintf('\nVan der Pol test 2\n')
          fn = 'f_VanderPol';
          Jn = 'J_VanderPol';
          fi = 'fi_VanderPol';
          fe = 'fe_VanderPol';
          Ji = 'Ji_VanderPol';
          Es = 'EStab_VanderPol';
          Tf = 12;
          tout = linspace(0,Tf,100);
          hmin = 1e-6;
          hmax = 1.0;
          atol = 1e-14*ones(2,1);
          global Pdata;
          Pdata.ep = 0.2;
          Pdata.ep2 = 0.2;
          u0 = 2;
          v0 = 0;
          Y0 = [u0; v0];

          % plot "true" solution
          if (show_plots == 1)
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
             [t,Ytrue] = ode15s('f_VanderPol', tout, Y0, opts);
             figure()
             plot(tout,Ytrue)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Van der Pol test 2','FontSize',14)
             set(gca,'FontSize',12)
             print('-depsc','vanderPol2')
          end

        
        case(6)  % van der Pol 3
          fprintf('\nVan der Pol test 3\n')
          fn = 'f_VanderPol';
          Jn = 'J_VanderPol';
          fi = 'fi_VanderPol';
          fe = 'fe_VanderPol';
          Ji = 'Ji_VanderPol';
          Es = 'EStab_VanderPol';
          Tf = 1000;
          tout = linspace(0,Tf,100);
          hmin = 1e-6;
          hmax = 10;
          atol = 1e-14*ones(2,1);
          global Pdata;
          Pdata.ep = 1e-3;
          Pdata.ep2 = 1e-3;
          u0 = 2;
          v0 = 0;
          Y0 = [u0; v0];
          
          % plot "true" solution
          if (show_plots == 1)
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
             [t,Ytrue] = ode15s('f_VanderPol', tout, Y0, opts);
             figure()
             plot(tout,Ytrue)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Van der Pol test 3','FontSize',14)
             set(gca,'FontSize',12)
             print('-depsc','vanderPol3')
          end
          
      end  % switch
      
      
      % loop over rtols
      for i_rtol = 1:length(rtols)

         % set relative tolerance for this sweep of tests
         rtol = rtols(i_rtol);
         fprintf('rtol = %g,  atol = %g\n',rtol,atol(1));

         
         % get "true" solution
         opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
         [t,Ytrue] = ode15s(fn, tout, Y0, opts);
         
         fprintf('\n');
         fprintf('method     ctrl   work   astep   sstep    maxerr   rmserr\n');
         fprintf('---------------------------------------------------------\n');
         
         % ode15s
         if (ode15s_meth)
            opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
            [t,Y] = ode15s(fn, [0,Tf], Y0, opts);  ns = length(t);
            [t,Y] = ode15s(fn, tout, Y0, opts);
            fprintf('ode15s       no    %i    %i    0    %.2e    %.2e\n',...
                    ns, ns, max(max(abs(Y-Ytrue))),...
                    sqrt(sum(sum((Y-Ytrue).^2))/numel(Y)));
         end

         % ode45
         if (ode45_meth)
            opts = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',hmin,'MaxStep',hmax);
            [t,Y] = ode45(fn, [0,Tf], Y0, opts);  ns = length(t);
            [t,Y] = ode45(fn, tout, Y0, opts);
            fprintf('ode45        no    %i    %i    0    %.2e    %.2e\n',...
                    ns, ns, max(max(abs(Y-Ytrue))),...
                    sqrt(sum(sum((Y-Ytrue).^2))/numel(Y)));
         end

         if (DIRK32_meth)
            mname = 'ARK3(2)4L[2]SA-ESDIRK';
            B = butcher(mname);

            if (I_cont)    % SDIRK3(2), I controller
               hmethod = 1;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK3(2)   I    %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (PI_cont)   % SDIRK3(2), PI controller
               hmethod = 2;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK3(2)   PI   %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (PID_cont)  % SDIRK3(2), PID controller
               hmethod = 3;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK3(2)   PID  %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end
            
            if (Gi_cont)   % SDIRK3(2), Gustafsson controller
               hmethod = 4;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK3(2)   iG   %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end
            
            if (Ge_cont)   % SDIRK3(2), explicit Gustafsson controller
               hmethod = 5;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK3(2)   eG   %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end
         end  % SDIRK32 method
         
         
         if (DIRK43_meth)
            mname = 'ARK4(3)6L[2]SA-ESDIRK';
            B = butcher(mname);
            
            if (I_cont)    % SDIRK4(3), I controller
               hmethod = 1;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK4(3)   I    %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end
            
            if (PI_cont)   % SDIRK4(3), PI controller
               hmethod = 2;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK4(3)   PI   %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
            end

            if (PID_cont)  % SDIRK4(3), PID controller
               hmethod = 3;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK4(3)   PID  %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs((Y'-Ytrue)./Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (Gi_cont)   % SDIRK4(3), Gustafsson controller
               hmethod = 4;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK4(3)   iG   %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (Ge_cont)   % SDIRK4(3), explicit Gustafsson controller
               hmethod = 5;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK4(3)   eG   %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end
         end  % SDIRK4(3) method
         
         
         if (DIRK54_meth)
            mname = 'ARK5(4)8L[2]SA-ESDIRK';
            B = butcher(mname);
            
            if (I_cont)    % SDIRK5(4), I controller
               hmethod = 1;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK5(4)   I    %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
            end

            if (PI_cont)   % SDIRK5(4), PI controller
               hmethod = 2;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK5(4)   PI   %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs((Y'-Ytrue)./Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (PID_cont)  % SDIRK5(4), PID controller
               hmethod = 3;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK5(4)   PID  %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (Gi_cont)   % SDIRK5(4), Gustafsson controller
               hmethod = 4;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK5(4)   iG   %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (Ge_cont)   % SDIRK5(4), explicit Gustafsson controller
               hmethod = 5;
               [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
               fprintf('SDIRK5(4)   eG   %i    %i    0    %.2e    %.2e\n',...
                       ns, ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end
         end  % SDIRK5(4) method


         if (ARK32_meth)
            mname = 'ARK3(2)4L[2]SA-ESDIRK';
            B = butcher(mname);
            mname2 = 'ARK3(2)4L[2]SA-ERK';
            B2 = butcher(mname2);

            if (I_cont)    % ARK3(2), I controller
               hmethod = 1;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK3(2)   I    %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (PI_cont)   % ARK3(2), PI controller
               hmethod = 2;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK3(2)   PI   %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (PID_cont)  % ARK3(2), PID controller
               hmethod = 3;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK3(2)   PID  %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (Gi_cont)   % ARK3(2), Gustafsson controller
               hmethod = 4;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK3(2)   iG   %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (Ge_cont)   % ARK3(2), explicit Gustafsson controller
               hmethod = 5;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK3(2)   eG   %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (Gie_cont)  % ARK3(2), imex Gustafsson controller
               hmethod = 6;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK3(2)   ieG  %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end
         end  % ARK3(2) method


         if (ARK43_meth)
            mname = 'ARK4(3)6L[2]SA-ESDIRK';
            B = butcher(mname);
            mname2 = 'ARK4(3)6L[2]SA-ERK';
            B2 = butcher(mname2);

            if (I_cont)    % ARK4(3), I controller
               hmethod = 1;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK4(3)   I   %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (PI_cont)   % ARK4(3), PI controller
               hmethod = 2;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK4(3)   PI  %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (PID_cont)  % ARK4(3), PID controller
               hmethod = 3;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK4(3)   PID %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (Gi_cont)   % ARK4(3), Gustafsson controller
               hmethod = 4;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK4(3)   iG   %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (Ge_cont)   % ARK4(3), explicit Gustafsson controller
               hmethod = 5;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK4(3)   eG   %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (Gie_cont)  % ARK4(3), imex Gustafsson controller
               hmethod = 6;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK4(3)   ieG  %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end
         end  % ARK4(3) method


         if (ARK54_meth)
            mname = 'ARK5(4)8L[2]SA-ESDIRK';
            B = butcher(mname);
            mname2 = 'ARK5(4)8L[2]SA-ERK';
            B2 = butcher(mname2);

            if (I_cont)   % ARK5(4), I controller
               hmethod = 1;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK5(4)   I    %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum(((Y'-Ytrue)./Ytrue).^2))/numel(Y)));
            end

            if (PI_cont)  % ARK5(4), PI controller
               hmethod = 2;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK5(4)   PI   %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs((Y'-Ytrue)./Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (PID_cont) % ARK5(4), PID controller
               hmethod = 3;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK5(4)   PID  %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (Gi_cont)  % ARK5(4), Gustafsson controller
               hmethod = 4;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK5(4)   iG   %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (Ge_cont)  % ARK5(4), explicit Gustafsson controller
               hmethod = 5;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK5(4)   eG   %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end

            if (Gie_cont) % ARK5(4), imex Gustafsson controller
               hmethod = 6;
               [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
               fprintf('ARK5(4)   ieG  %i    0    0    %.2e    %.2e\n',...
                       ns, max(max(abs(Y'-Ytrue))),...
                       sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
            end
         end  % ARK5(4) method

         fprintf('---------------------------------------------------------\n');
         fprintf('\n\n')
         
      end  % i_rtol loop

   end  % if tests(i_test)
   
end  % i_test loop
   

% stop diary
disp('   ')
diary off

% exit Matlab
exit

% end of testing script