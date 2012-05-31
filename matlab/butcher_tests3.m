% driver to check accuracy of a few specific RK methods (and embeddings) in butcher.m
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% May 2012
% All Rights Reserved

clear
!\rm butcher_tests3.txt
diary butcher_tests3.txt


% ARK3(2)4L[2]SA-ERK method
mname = 'ARK3(2)4L[2]SA-ERK';
B = butcher(mname);
butcher_test_method(B,mname);

% ARK3(2)4L[2]SA-ESDIRK method
mname = 'ARK3(2)4L[2]SA-ESDIRK';
B = butcher(mname);
butcher_test_method(B,mname);

% ARK4(3)6L[2]SA-ERK method
mname = 'ARK4(3)6L[2]SA-ERK';
B = butcher(mname);
butcher_test_method(B,mname);

% ARK4(3)6L[2]SA-ESDIRK method
mname = 'ARK4(3)6L[2]SA-ESDIRK';
B = butcher(mname);
butcher_test_method(B,mname);

% ARK5(4)8L[2]SA-ERK method
mname = 'ARK5(4)8L[2]SA-ERK';
B = butcher(mname);
butcher_test_method(B,mname);

% ARK5(4)8L[2]SA-ESDIRK method
mname = 'ARK5(4)8L[2]SA-ESDIRK';
B = butcher(mname);
butcher_test_method(B,mname);

% Sayfy-Aburub-4-3-ERK method
mname = 'Sayfy-Aburub-4-3-ERK';
B = butcher(mname);
butcher_test_method(B,mname);

% Sayfy-Aburub-4-3-DIRK method
mname = 'Sayfy-Aburub-4-3-DIRK';
B = butcher(mname);
butcher_test_method(B,mname);

% Heun-Euler-ERK method
mname = 'Heun-Euler-ERK';
B = butcher(mname);
butcher_test_method(B,mname);

% Bogacki-Shampine-ERK method
mname = 'Bogacki-Shampine-ERK';
B = butcher(mname);
butcher_test_method(B,mname);

% Fehlberg-ERK method
mname = 'Fehlberg-ERK';
B = butcher(mname);
butcher_test_method(B,mname);

% Cash-Karp-ERK method
mname = 'Cash-Karp-ERK';
B = butcher(mname);
butcher_test_method(B,mname);

% Dormand-Prince-ERK method
mname = 'Dormand-Prince-ERK';
B = butcher(mname);
butcher_test_method(B,mname);

% Merson-4-5 method
mname = 'Merson-4-5-ERK';
B = butcher(mname);
butcher_test_method(B,mname);

% Zonneveld-4-3-ERK method
mname = 'Zonneveld-4-3-ERK';
B = butcher(mname);
butcher_test_method(B,mname);

% Verner-6-5-ERK method
mname = 'Verner-6-5-ERK';
B = butcher(mname);
butcher_test_method(B,mname);

% Fehlberg-8-7-ERK method
mname = 'Fehlberg-8-7-ERK';
B = butcher(mname);
butcher_test_method(B,mname);

% TRBDF2-ESDIRK method
mname = 'TRBDF2-ESDIRK';
B = butcher(mname);
butcher_test_method(B,mname);

% TRX2-ESDIRK method
mname = 'TRX2-ESDIRK';
B = butcher(mname);
butcher_test_method(B,mname);

% Billington-SDIRK method
mname = 'Billington-SDIRK';
B = butcher(mname);
butcher_test_method(B,mname);

% Cash(5,2,4)-SDIRK method
mname = 'Cash(5,2,4)-SDIRK';
B = butcher(mname);
butcher_test_method(B,mname);

% Cash(5,3,4)-SDIRK method
mname = 'Cash(5,3,4)-SDIRK';
B = butcher(mname);
butcher_test_method(B,mname);

% SDIRK-4-5 method
mname = 'SDIRK-4-5';
B = butcher(mname);
butcher_test_method(B,mname);

% Kvaerno(4,2,3)-ESDIRK method
mname = 'Kvaerno(4,2,3)-ESDIRK';
B = butcher(mname);
butcher_test_method(B,mname);

% Kvaerno(5,3,4)-ESDIRK method
mname = 'Kvaerno(5,3,4)-ESDIRK';
B = butcher(mname);
butcher_test_method(B,mname);

% Kvaerno(7,4,5)-ESDIRK method
mname = 'Kvaerno(7,4,5)-ESDIRK';
B = butcher(mname);
butcher_test_method(B,mname);

% Ismail(7,4,5)-ESDIRK method
mname = 'Ismail(7,4,5)-ESDIRK';
B = butcher(mname);
butcher_test_method(B,mname);

diary off
% end of script
