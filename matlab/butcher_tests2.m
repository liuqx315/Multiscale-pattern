% driver to check analytical properties of various RK methods in butcher.m
%
% Daniel R. Reynolds
% SMU Mathematics
% May 2012
% All Rights Reserved

clear
!\rm butcher_tests2.txt
diary butcher_tests2.txt


fprintf('          Method          q   p   Bs   As   Ls\n');
fprintf('-----------------------------------------------\n');


%%%%%%%%%%%%%%%%%%% ERK Methods %%%%%%%%%%%%%%%%%%%

%    ARK3(2)4L[2]SA-ERK
mname = 'ARK3(2)4L[2]SA-ERK';
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

%    ARK5(4)8L[2]SA-ERK
mname = 'ARK5(4)8L[2]SA-ERK';
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

%    Ascher(2,3,3)-ERK
mname = 'Ascher(2,3,3)-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Ascher(2,3,2)-ERK
mname = 'Ascher(2,3,2)-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Ascher(2,2,2)-ERK
mname = 'Ascher(2,2,2)-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Ascher(3,4,3)-ERK
mname = 'Ascher(3,4,3)-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Ascher(4,4,3)-ERK
mname = 'Ascher(4,4,3)-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Cooper4-ERK
mname = 'Cooper4-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Cooper6-ERK
mname = 'Cooper6-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

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

%    Fehlberg-ERK
mname = 'Fehlberg-ERK';
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

%    Dormand-Prince-ERK
mname = 'Dormand-Prince-ERK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    ERK-1-1
mname = 'ERK-1-1';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    ERK-2-2
mname = 'ERK-2-2';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    ERK-3-3
mname = 'ERK-3-3';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    ERK-4-4
mname = 'ERK-4-4';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Merson-4-5-ERK
mname = 'Merson-4-5-ERK';
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


%%%%%%%%%%%%%%%%%%% DIRK Methods %%%%%%%%%%%%%%%%%%%


%    ARK3(2)4L[2]SA-ESDIRK
mname = 'ARK3(2)4L[2]SA-ESDIRK';
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

%    ARK5(4)8L[2]SA-ESDIRK
mname = 'ARK5(4)8L[2]SA-ESDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Sayfy-Aburub-4-3-DIRK
mname = 'Sayfy-Aburub-4-3-DIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Ascher(2,3,3)-SDIRK
mname = 'Ascher(2,3,3)-SDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Ascher(3,4,3)-SDIRK
mname = 'Ascher(3,4,3)-SDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Ascher(2,3,2)-SDIRK
mname = 'Ascher(2,3,2)-SDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Ascher(2,2,2)-SDIRK
mname = 'Ascher(2,2,2)-SDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Ascher(4,4,3)-SDIRK
mname = 'Ascher(4,4,3)-SDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Cooper4-ESDIRK
mname = 'Cooper4-ESDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Cooper6-ESDIRK
mname = 'Cooper6-ESDIRK';
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

%    TRX2-ESDIRK
mname = 'TRX2-ESDIRK';
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

%    Ismail(7,4,5)-ESDIRK
mname = 'Ismail(7,4,5)-ESDIRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    SDIRK-2-2
mname = 'SDIRK-2-2';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    SDIRK-4-5
mname = 'SDIRK-4-5';
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

%    Kvaerno(5,3,4)-ESDIRK 
mname = 'Kvaerno(5,3,4)-ESDIRK';
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



%%%%%%%%%%%%%%%%%%% IRK Methods %%%%%%%%%%%%%%%%%%%

%    IRK-1-1
mname = 'IRK-1-1';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Crank-Nicolson-2-2-IRK
mname = 'Crank-Nicolson-2-2-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    SIRK-2-2
mname = 'SIRK-2-2';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    ESIRK-2-2
mname = 'ESIRK-2-2';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Gauss-2-4-IRK
mname = 'Gauss-2-4-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    RadauIIA-2-3-IRK
mname = 'RadauIIA-2-3-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIII-2-2-IRK
mname = 'LobattoIII-2-2-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIIIA-2-2-IRK
mname = 'LobattoIIIA-2-2-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIIIB-2-2-IRK
mname = 'LobattoIIIB-2-2-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIIIC-2-2-IRK
mname = 'LobattoIIIC-2-2-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Gauss-3-2-IRK
mname = 'Gauss-3-2-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    RadauI-3-5-IRK
mname = 'RadauI-3-5-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    RadauIA-3-5-IRK
mname = 'RadauIA-3-5-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    RadauII-3-5-IRK
mname = 'RadauII-3-5-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    RadauIIA-3-5-IRK
mname = 'RadauIIA-3-5-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIII-3-4-IRK
mname = 'LobattoIII-3-4-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIIIA-3-4-IRK
mname = 'LobattoIIIA-3-4-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIIIB-3-4-IRK
mname = 'LobattoIIIB-3-4-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIIIC-3-4-IRK
mname = 'LobattoIIIC-3-4-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIII-4-6-IRK
mname = 'LobattoIII-4-6-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIIIA-4-6-IRK
mname = 'LobattoIIIA-4-6-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIIIB-4-6-IRK
mname = 'LobattoIIIB-4-6-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIIIC-4-6-IRK
mname = 'LobattoIIIC-4-6-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    RadauIIA-5-9-IRK
mname = 'RadauIIA-5-9-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIII-5-8-IRK
mname = 'LobattoIII-5-8-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIIIA-5-8-IRK
mname = 'LobattoIIIA-5-8-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIIIB-5-8-IRK
mname = 'LobattoIIIB-5-8-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    LobattoIIIC-5-8-IRK
mname = 'LobattoIIIC-5-8-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Gauss-6-12-IRK
mname = 'Gauss-6-12-IRK';
fprintf(' %22s',mname);
B = butcher(mname);
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('  %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

%    Qualifying exam
mname = 'Qualifying Exam';
fprintf(' %22s',mname);
B = [0, 0, 0, 0; 1/2, 1/4, 1/4, 0; 1, 1/4, 1/2, 1/4; 0, 1/6, 2/3, 1/6];
[q,p,Bs,As,Ls] = check_butcher(B);
fprintf('   %i   %i   %i    %i    %i\n',q,p,Bs,As,Ls);

fprintf('-----------------------------------------------\n');

diary off
% end of script
