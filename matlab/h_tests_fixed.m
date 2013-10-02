%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------
% driver that runs full suite of test problems to examine different time
% step selection strategies
clear
maxNumCompThreads(1);

% tests to run
b1_test = 0;  % brusselator test 1
b2_test = 0;  % brusselator test 2
b3_test = 0;  % brusselator test 3
v1_test = 1;  % van der Pol test 1
v2_test = 1;  % van der Pol test 2
v3_test = 1;  % van der Pol test 3
a1_test = 1;  % analytical  test 1
a2_test = 1;  % analytical  test 2
a3_test = 1;  % analytical  test 3
tests = [b1_test, b2_test, b3_test, v1_test, v2_test, v3_test, a1_test, a2_test, a3_test];

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

% set controller to PID (unused since this test uses fixed steps)
hmethod = 3;

% flag to generate plots
show_plots = 1;

% relative step sizes to use
hfacs = [1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4];

%    start diary
!\rm -f h_test_results_fixed2.txt
diary h_test_results_fixed2.txt


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
          tout = linspace(0,Tf,101);
          rtol = 1e-2;
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
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',1e-8, 'MaxStep',1e-2);
             [t,Ytrue] = ode15s('f_bruss', tout, Y0, opts);
             figure()
             plot(tout,Ytrue,'LineWidth',2)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Brusselator test 1','FontSize',14)
             legend('u','v','w')
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
          tout = linspace(0,Tf,101);
          rtol = 1e-2;
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
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',1e-8, 'MaxStep',1e-2);
             [t,Ytrue] = ode15s('f_bruss', tout, Y0, opts);
             figure()
             plot(tout,Ytrue,'LineWidth',2)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Brusselator test 2','FontSize',14)
             legend('u','v','w')
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
          tout = linspace(0,Tf,101);
          rtol = 1e-2;
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
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',1e-8, 'MaxStep',1e-2);
             [t,Ytrue] = ode15s('f_bruss', tout, Y0, opts);
             figure()
             plot(tout,Ytrue,'LineWidth',2)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Brusselator test 3','FontSize',14)
             legend('u','v','w')
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
          tout = linspace(0,Tf,101);
          rtol = 1e-2;
          atol = 1e-14*ones(2,1);
          global Pdata;
          Pdata.ep = 1e-3;
          Pdata.ep2 = 1.0;
          u0 = -2;
          v0 = -2.3553013976081;
          Y0 = [u0; v0];

          % plot "true" solution
          if (show_plots == 1)
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',1e-8, 'MaxStep',1e-2);
             [t,Ytrue] = ode15s('f_VanderPol', tout, Y0, opts);
             figure()
             plot(tout,Ytrue,'LineWidth',2)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Van der Pol test 1','FontSize',14)
             legend('u','v')
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
          tout = linspace(0,Tf,101);
          rtol = 1e-2;
          atol = 1e-14*ones(2,1);
          global Pdata;
          Pdata.ep = 0.2;
          Pdata.ep2 = 0.2;
          u0 = 2;
          v0 = 0;
          Y0 = [u0; v0];

          % plot "true" solution
          if (show_plots == 1)
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',1e-8, 'MaxStep',1e-2);
             [t,Ytrue] = ode15s('f_VanderPol', tout, Y0, opts);
             figure()
             plot(tout,Ytrue,'LineWidth',2)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Van der Pol test 2','FontSize',14)
             legend('u','v',2)
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
          %Tf = 1000;
          Tf = 800;
          tout = linspace(0,Tf,101);
          rtol = 1e-2;
          atol = 1e-14*ones(2,1);
          global Pdata;
          Pdata.ep = 1e-3;
          Pdata.ep2 = 1e-3;
          u0 = 2;
          v0 = 0;
          Y0 = [u0; v0];
          
          % plot "true" solution
          if (show_plots == 1)
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',1e-8, 'MaxStep',1e-2);
             [t,Ytrue] = ode15s('f_VanderPol', tout, Y0, opts);
             figure()
             plot(tout,Ytrue,'LineWidth',2)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Van der Pol test 3','FontSize',14)
             legend('u','v')
             set(gca,'FontSize',12)
             print('-depsc','vanderPol3')
          end

        
        case(7)  % analytical test 1
          fprintf('\nAnalytical\n')
          fn = 'f_analytic';
          Jn = 'J_analytic';
          fi = 'fi_analytic';
          fe = 'fe_analytic';
          Ji = 'Ji_analytic';
          Es = 'EStab_analytic';
          Tf = 10;
          tout = linspace(0,Tf,101);
          rtol = 1e-2;
          atol = 1e-14*ones(1,1);
          global Pdata;
          Pdata.ep = 1e-0;
          Y0 = 2;
          
          % plot "true" solution
          if (show_plots == 1)
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',1e-8, 'MaxStep',1e-2);
             [t,Ytrue] = ode15s('f_analytic', tout, Y0, opts);
             figure()
             plot(tout,Ytrue,'LineWidth',2)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Variable stiffness Test 1','FontSize',14)
             axis([0 10 0 2])
             set(gca,'FontSize',12)
             print('-depsc','analytic1')
          end

        
        case(8)  % analytical test 2
          fprintf('\nAnalytical\n')
          fn = 'f_analytic';
          Jn = 'J_analytic';
          fi = 'fi_analytic';
          fe = 'fe_analytic';
          Ji = 'Ji_analytic';
          Es = 'EStab_analytic';
          Tf = 10;
          tout = linspace(0,Tf,101);
          rtol = 1e-2;
          atol = 1e-14*ones(1,1);
          global Pdata;
          Pdata.ep = 1e-1;
          Y0 = 2;
          
          % plot "true" solution
          if (show_plots == 1)
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',1e-8, 'MaxStep',1e-2);
             [t,Ytrue] = ode15s('f_analytic', tout, Y0, opts);
             figure()
             plot(tout,Ytrue,'LineWidth',2)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Variable stiffness Test 2','FontSize',14)
             axis([0 10 0 2])
             set(gca,'FontSize',12)
             print('-depsc','analytic2')
          end

        
        case(9)  % analytical test 3
          fprintf('\nAnalytical\n')
          fn = 'f_analytic';
          Jn = 'J_analytic';
          fi = 'fi_analytic';
          fe = 'fe_analytic';
          Ji = 'Ji_analytic';
          Es = 'EStab_analytic';
          Tf = 10;
          tout = linspace(0,Tf,101);
          rtol = 1e-2;
          atol = 1e-14*ones(1,1);
          global Pdata;
          Pdata.ep = 1e-2;
          Y0 = 2;
          
          % plot "true" solution
          if (show_plots == 1)
             opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',1e-8, 'MaxStep',1e-2);
             [t,Ytrue] = ode15s('f_analytic', tout, Y0, opts);
             figure()
             plot(tout,Ytrue,'LineWidth',2)
             xlabel('t','FontSize',12), ylabel('y','FontSize',12)
             title('Variable stiffness Test 3','FontSize',14)
             axis([0 10 0 2])
             set(gca,'FontSize',12)
             print('-depsc','analytic3')
          end
          
      end  % switch
      
      
      % loop over hfacs
      for i_hfac = 1:length(hfacs)

         % set time step for this sweep of tests
         hfac = hfacs(i_hfac);
         fprintf('h_factor = %g\n',hfac);
         hmin = hfac*Tf;
         hmax = hfac*Tf;
         
         % get "true" solution
         opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',1e-8, 'MaxStep',hmax);
         [t,Ytrue] = ode15s(fn, tout, Y0, opts);
         
         fprintf('\n');
         fprintf('method     h    maxerr   rmserr\n');
         fprintf('-------------------------------\n');
         
         % ode15s
         if (ode15s_meth)
            opts = odeset('InitialStep',hmin,'MaxStep',hmax);
            [t,Y] = ode15s(fn, [0,Tf], Y0, opts);  ns = length(t);
            [t,Y] = ode15s(fn, tout, Y0, opts);
            fprintf('ode15s   %g    %.2e    %.2e\n', hmax, ...
                    max(max(abs(Y-Ytrue))), sqrt(sum(sum((Y-Ytrue).^2))/numel(Y)));
         end

         % ode45
         if (ode45_meth)
            opts = odeset('InitialStep',hmin,'MaxStep',hmax);
            [t,Y] = ode45(fn, [0,Tf], Y0, opts);  ns = length(t);
            [t,Y] = ode45(fn, tout, Y0, opts);
            fprintf('ode45    %g    %.2e    %.2e\n', hmax, ...
                    max(max(abs(Y-Ytrue))), sqrt(sum(sum((Y-Ytrue).^2))/numel(Y)));
         end

         % SDIRK3(2) 
         if (DIRK32_meth)
            mname = 'ARK3(2)4L[2]SA-ESDIRK';
            B = butcher(mname);
            [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
            fprintf('SDIRK3(2)   %g    %.2e    %.2e\n', hmax, ...
                    max(max(abs(Y'-Ytrue))), sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
         end
         
         % SDIRK4(3)
         if (DIRK43_meth)
            mname = 'ARK4(3)6L[2]SA-ESDIRK';
            B = butcher(mname);
            [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
            fprintf('SDIRK4(3)   %g    %.2e    %.2e\n', hmax, ...
                    max(max(abs(Y'-Ytrue))), sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
         end
         
         % SDIRK5(4)
         if (DIRK54_meth)
            mname = 'ARK5(4)8L[2]SA-ESDIRK';
            B = butcher(mname);
            [t,Y,ns] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmethod);
            fprintf('SDIRK5(4)   %g    %.2e    %.2e\n', hmax, ...
                    max(max(abs(Y'-Ytrue))), sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
         end

         % ARK3(2)
         if (ARK32_meth)
            mname = 'ARK3(2)4L[2]SA-ESDIRK';
            B = butcher(mname);
            mname2 = 'ARK3(2)4L[2]SA-ERK';
            B2 = butcher(mname2);
            [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
            fprintf('ARK3(2)   %g    %.2e    %.2e\n', hmax, ...
                    max(max(abs(Y'-Ytrue))), sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
         end

         % ARK4(3)
         if (ARK43_meth)
            mname = 'ARK4(3)6L[2]SA-ESDIRK';
            B = butcher(mname);
            mname2 = 'ARK4(3)6L[2]SA-ERK';
            B2 = butcher(mname2);
            [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
            fprintf('ARK4(3)   %g    %.2e    %.2e\n', hmax, ...
                    max(max(abs(Y'-Ytrue))), sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
         end

         % ARK5(4)
         if (ARK54_meth)
            mname = 'ARK5(4)8L[2]SA-ESDIRK';
            B = butcher(mname);
            mname2 = 'ARK5(4)8L[2]SA-ERK';
            B2 = butcher(mname2);
            [t,Y,ns] = solve_ARK(fi, fe, Ji, Es, tout, Y0, B, B2, rtol, atol, hmin, hmax, hmethod);
            fprintf('ARK5(4)   %g    %.2e    %.2e\n', hmax, ...
                    max(max(abs(Y'-Ytrue))), sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y)));
         end

         fprintf('---------------------------------------------------------\n');
         fprintf('\n\n')
         
      end  % i_hfac loop

   end  % if tests(i_test)
   
end  % i_test loop
   

% stop diary
disp('   ')
diary off

% exit Matlab
%exit

% end of testing script