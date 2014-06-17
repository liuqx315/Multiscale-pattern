% driver to check accuracy of a few specific RK methods (and embeddings) in butcher.m
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

clear
!\rm verify_butcher.txt
diary verify_butcher.txt


% first, examine the butcher table itself
disp('  ')
disp('Analyzing Butcher tables themselves:')
disp('  ')


fprintf('          Method          q   p   Bs   As   Ls\n');
fprintf('-----------------------------------------------\n');

%    Heun-Euler-ERK
mname = 'Heun-Euler-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Bogacki-Shampine-ERK
mname = 'Bogacki-Shampine-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    ARK3(2)4L[2]SA-ERK
mname = 'ARK3(2)4L[2]SA-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Zonneveld-4-3-ERK
mname = 'Zonneveld-4-3-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    ARK4(3)6L[2]SA-ERK
mname = 'ARK4(3)6L[2]SA-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Sayfy-Aburub-4-3-ERK
mname = 'Sayfy-Aburub-4-3-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Cash-Karp-ERK
mname = 'Cash-Karp-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Fehlberg-ERK
mname = 'Fehlberg-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Dormand-Prince-ERK
mname = 'Dormand-Prince-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    ARK5(4)8L[2]SA-ERK
mname = 'ARK5(4)8L[2]SA-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Verner-6-5-ERK
mname = 'Verner-6-5-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Fehlberg-8-7-ERK
mname = 'Fehlberg-8-7-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);




%    SDIRK-2-1
mname = 'SDIRK-2-1';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Billington-SDIRK
mname = 'Billington-SDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    TRBDF2-ESDIRK
mname = 'TRBDF2-ESDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Kvaerno(4,2,3)-ESDIRK 
mname = 'Kvaerno(4,2,3)-ESDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    ARK3(2)4L[2]SA-ESDIRK
mname = 'ARK3(2)4L[2]SA-ESDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Cash(5,2,4)-SDIRK
mname = 'Cash(5,2,4)-SDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Cash(5,3,4)-SDIRK
mname = 'Cash(5,3,4)-SDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    SDIRK-5-4
mname = 'SDIRK-5-4';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Kvaerno(5,3,4)-ESDIRK 
mname = 'Kvaerno(5,3,4)-ESDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    ARK4(3)6L[2]SA-ESDIRK
mname = 'ARK4(3)6L[2]SA-ESDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Kvaerno(7,4,5)-ESDIRK 
mname = 'Kvaerno(7,4,5)-ESDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    ARK5(4)8L[2]SA-ESDIRK
mname = 'ARK5(4)8L[2]SA-ESDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

fprintf('-----------------------------------------------\n\n');





%------------------------------------------
% second, estimate the numerical accuracy of each method on a
% simple test problem
disp('  ')
disp('Estimating the numerical accuracy of the methods,')
disp('  first is the main method, second the embedding:')
disp('  ')


% Heun-Euler-ERK method
mname = 'Heun-Euler-ERK';
B = butcher(mname);  B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Heun-Euler-ERK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Bogacki-Shampine-ERK method
mname = 'Bogacki-Shampine-ERK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Bogacki-Shampine-ERK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% ARK3(2)4L[2]SA-ERK method
mname = 'ARK3(2)4L[2]SA-ERK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'ARK3(2)4L[2]SA-ERK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Zonneveld-4-3-ERK method
mname = 'Zonneveld-4-3-ERK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Zonneveld-4-3-ERK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% ARK4(3)6L[2]SA-ERK method
mname = 'ARK4(3)6L[2]SA-ERK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'ARK4(3)6L[2]SA-ERK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Sayfy-Aburub-4-3-ERK method
mname = 'Sayfy-Aburub-4-3-ERK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Sayfy-Aburub-4-3-ERK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Cash-Karp-ERK method
mname = 'Cash-Karp-ERK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Cash-Karp-ERK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Fehlberg-ERK method
mname = 'Fehlberg-ERK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Fehlberg-ERK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Dormand-Prince-ERK method
mname = 'Dormand-Prince-ERK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Dormand-Prince-ERK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% ARK5(4)8L[2]SA-ERK method
mname = 'ARK5(4)8L[2]SA-ERK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'ARK5(4)8L[2]SA-ERK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Verner-6-5-ERK method
mname = 'Verner-6-5-ERK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Verner-6-5-ERK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Fehlberg-8-7-ERK method
mname = 'Fehlberg-8-7-ERK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Fehlberg-8-7-ERK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);



% SDIRK-2-1 method
mname = 'SDIRK-2-1';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'SDIRK-2-1';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Billington-SDIRK method
mname = 'Billington-SDIRK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Billington-SDIRK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% TRBDF2-ESDIRK method
mname = 'TRBDF2-ESDIRK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'TRBDF2-ESDIRK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Kvaerno(4,2,3)-ESDIRK method
mname = 'Kvaerno(4,2,3)-ESDIRK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Kvaerno(4,2,3)-ESDIRK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% ARK3(2)4L[2]SA-ESDIRK method
mname = 'ARK3(2)4L[2]SA-ESDIRK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'ARK3(2)4L[2]SA-ESDIRK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Cash(5,2,4)-SDIRK method
mname = 'Cash(5,2,4)-SDIRK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Cash(5,2,4)-SDIRK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Cash(5,3,4)-SDIRK method
mname = 'Cash(5,3,4)-SDIRK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Cash(5,3,4)-SDIRK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% SDIRK-5-4 method
mname = 'SDIRK-5-4';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'SDIRK-5-4';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Kvaerno(5,3,4)-ESDIRK method
mname = 'Kvaerno(5,3,4)-ESDIRK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Kvaerno(5,3,4)-ESDIRK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% ARK4(3)6L[2]SA-ESDIRK method
mname = 'ARK4(3)6L[2]SA-ESDIRK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'ARK4(3)6L[2]SA-ESDIRK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% Kvaerno(7,4,5)-ESDIRK method
mname = 'Kvaerno(7,4,5)-ESDIRK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'Kvaerno(7,4,5)-ESDIRK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);

% ARK5(4)8L[2]SA-ESDIRK method
mname = 'ARK5(4)8L[2]SA-ESDIRK';
B = butcher(mname); B = B(1:end-1,:);
butcher_test_method(B,mname);
mname = 'ARK5(4)8L[2]SA-ESDIRK';
B = butcher(mname);  B = [B(1:end-2,:); B(end,:)];
butcher_test_method(B,mname);



diary off
% end of script
