% script file that goes through all methods in butcher.m and builds
% coefficients for dense output (uses new approach for
% coefficients)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved
clear
format long

% start output diary
!\rm -f dense_tables.txt
diary dense_tables.txt

% go through methods, outputting dense output matrix for each to screen
B = build_dense('ARK3(2)4L[2]SA-ERK');
fprintf('\n\nARK3(2)4L[2]SA-ERK dense output coefficients:\n')
disp(B)

B = build_dense('ARK3(2)4L[2]SA-ESDIRK');
fprintf('\n\nARK3(2)4L[2]SA-ESDIRK dense output coefficients:\n')
disp(B)

B = build_dense('ARK4(3)6L[2]SA-ERK');
fprintf('\n\nARK4(3)6L[2]SA-ERK dense output coefficients:\n')
disp(B)

B = build_dense('ARK4(3)6L[2]SA-ESDIRK');
fprintf('\n\nARK4(3)6L[2]SA-ESDIRK dense output coefficients:\n')
disp(B)

B = build_dense('ARK5(4)8L[2]SA-ERK');
fprintf('\n\nARK5(4)8L[2]SA-ERK dense output coefficients:\n')
disp(B)

B = build_dense('ARK5(4)8L[2]SA-ESDIRK');
fprintf('\n\nARK5(4)8L[2]SA-ESDIRK dense output coefficients:\n')
disp(B)

B = build_dense('Ascher(2,3,3)-ERK');
fprintf('\n\nAscher(2,3,3)-ERK dense output coefficients:\n')
disp(B)

B = build_dense('Ascher(2,3,3)-SDIRK');
fprintf('\n\nAscher(2,3,3)-SDIRK dense output coefficients:\n')
disp(B)

B = build_dense('Ascher(2,3,2)-ERK');
fprintf('\n\nAscher(2,3,2)-ERK dense output coefficients:\n')
disp(B)

B = build_dense('Ascher(2,3,2)-SDIRK');
fprintf('\n\nAscher(2,3,2)-SDIRK dense output coefficients:\n')
disp(B)

B = build_dense('Ascher(2,2,2)-ERK');
fprintf('\n\nAscher(2,2,2)-ERK dense output coefficients:\n')
disp(B)

B = build_dense('Ascher(2,2,2)-SDIRK');
fprintf('\n\nAscher(2,2,2)-SDIRK dense output coefficients:\n')
disp(B)

B = build_dense('Ascher(3,4,3)-ERK');
fprintf('\n\nAscher(3,4,3)-ERK dense output coefficients:\n')
disp(B)

B = build_dense('Ascher(3,4,3)-SDIRK');
fprintf('\n\nAscher(3,4,3)-SDIRK dense output coefficients:\n')
disp(B)

B = build_dense('Ascher(4,4,3)-ERK');
fprintf('\n\nAscher(4,4,3)-ERK dense output coefficients:\n')
disp(B)

B = build_dense('Ascher(4,4,3)-SDIRK');
fprintf('\n\nAscher(4,4,3)-SDIRK dense output coefficients:\n')
disp(B)

B = build_dense('Cooper4-ERK');
fprintf('\n\nCooper4-ERK dense output coefficients:\n')
disp(B)

B = build_dense('Cooper4-ESDIRK');
fprintf('\n\nCooper4-ESDIRK dense output coefficients:\n')
disp(B)

B = build_dense('Cooper6-ERK');
fprintf('\n\nCooper6-ERK dense output coefficients:\n')
disp(B)

B = build_dense('Cooper6-ESDIRK');
fprintf('\n\nCooper6-ESDIRK dense output coefficients:\n')
disp(B)

B = build_dense('Heun-Euler-ERK');
fprintf('\n\nHeun-Euler-ERK dense output coefficients:\n')
disp(B)

B = build_dense('Bogacki-Shampine-ERK');
fprintf('\n\nBogacki-Shampine-ERK dense output coefficients:\n')
disp(B)

B = build_dense('Fehlberg-ERK');
fprintf('\n\nFehlberg-ERK dense output coefficients:\n')
disp(B)

B = build_dense('Cash-Karp-ERK');
fprintf('\n\nCash-Karp-ERK dense output coefficients:\n')
disp(B)

B = build_dense('Dormand-Prince-ERK');
fprintf('\n\nDormand-Prince-ERK dense output coefficients:\n')
disp(B)

B = build_dense('ERK-1-1');
fprintf('\n\nERK-1-1 dense output coefficients:\n')
disp(B)

B = build_dense('ERK-2-2');
fprintf('\n\nERK-2-2 dense output coefficients:\n')
disp(B)

B = build_dense('ERK-3-3');
fprintf('\n\nERK-3-3 dense output coefficients:\n')
disp(B)

B = build_dense('ERK-4-4');
fprintf('\n\nERK-4-4 dense output coefficients:\n')
disp(B)

B = build_dense('TRBDF2-ESDIRK');
fprintf('\n\nTRBDF2-ESDIRK dense output coefficients:\n')
disp(B)

B = build_dense('TRX2-ESDIRK');
fprintf('\n\nTRX2-ESDIRK dense output coefficients:\n')
disp(B)

B = build_dense('Billington-SDIRK');
fprintf('\n\nBillington-SDIRK dense output coefficients:\n')
disp(B)

B = build_dense('Cash(5,2,4)-SDIRK');
fprintf('\n\nCash(5,2,4)-SDIRK dense output coefficients:\n')
disp(B)

B = build_dense('Cash(5,3,4)-SDIRK');
fprintf('\n\nCash(5,3,4)-SDIRK dense output coefficients:\n')
disp(B)

B = build_dense('Ismail(7,4,5)-ESDIRK');
fprintf('\n\nIsmail(7,4,5)-ESDIRK dense output coefficients:\n')
disp(B)

B = build_dense('SDIRK-2-2');
fprintf('\n\nSDIRK-2-2 dense output coefficients:\n')
disp(B)

B = build_dense('IRK-1-1');
fprintf('\n\nIRK-1-1 dense output coefficients:\n')
disp(B)

B = build_dense('Crank-Nicolson-2-2-IRK');
fprintf('\n\nCrank-Nicolson-2-2-IRK dense output coefficients:\n')
disp(B)

B = build_dense('SIRK-2-2');
fprintf('\n\nSIRK-2-2 dense output coefficients:\n')
disp(B)

B = build_dense('ESIRK-2-2');
fprintf('\n\nESIRK-2-2 dense output coefficients:\n')
disp(B)

B = build_dense('Gauss-2-4-IRK');
fprintf('\n\nGauss-2-4-IRK dense output coefficients:\n')
disp(B)

B = build_dense('RadauIIA-2-3-IRK');
fprintf('\n\nRadauIIA-2-3-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIII-2-2-IRK');
fprintf('\n\nLobattoIII-2-2-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIIIA-2-2-IRK');
fprintf('\n\nLobattoIIIA-2-2-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIIIB-2-2-IRK');
fprintf('\n\nLobattoIIIB-2-2-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIIIC-2-2-IRK');
fprintf('\n\nLobattoIIIC-2-2-IRK dense output coefficients:\n')
disp(B)

B = build_dense('Gauss-3-2-IRK');
fprintf('\n\nGauss-3-2-IRK dense output coefficients:\n')
disp(B)

B = build_dense('RadauI-3-5-IRK');
fprintf('\n\nRadauI-3-5-IRK dense output coefficients:\n')
disp(B)

B = build_dense('RadauIA-3-5-IRK');
fprintf('\n\nRadauIA-3-5-IRK dense output coefficients:\n')
disp(B)

B = build_dense('RadauII-3-5-IRK');
fprintf('\n\nRadauII-3-5-IRK dense output coefficients:\n')
disp(B)

B = build_dense('RadauIIA-3-5-IRK');
fprintf('\n\nRadauIIA-3-5-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIII-3-4-IRK');
fprintf('\n\nLobattoIII-3-4-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIIIA-3-4-IRK');
fprintf('\n\nLobattoIIIA-3-4-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIIIB-3-4-IRK');
fprintf('\n\nLobattoIIIB-3-4-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIIIC-3-4-IRK');
fprintf('\n\nLobattoIIIC-3-4-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIII-4-6-IRK');
fprintf('\n\nLobattoIII-4-6-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIIIA-4-6-IRK');
fprintf('\n\nLobattoIIIA-4-6-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIIIB-4-6-IRK');
fprintf('\n\nLobattoIIIB-4-6-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIIIC-4-6-IRK');
fprintf('\n\nLobattoIIIC-4-6-IRK dense output coefficients:\n')
disp(B)

B = build_dense('RadauIIA-5-9-IRK');
fprintf('\n\nRadauIIA-5-9-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIII-5-8-IRK');
fprintf('\n\nLobattoIII-5-8-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIIIA-5-8-IRK');
fprintf('\n\nLobattoIIIA-5-8-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIIIB-5-8-IRK');
fprintf('\n\nLobattoIIIB-5-8-IRK dense output coefficients:\n')
disp(B)

B = build_dense('LobattoIIIC-5-8-IRK');
fprintf('\n\nLobattoIIIC-5-8-IRK dense output coefficients:\n')
disp(B)

B = build_dense('SDIRK-4-5');
fprintf('\n\nSDIRK-4-5 dense output coefficients:\n')
disp(B)

B = build_dense('Gauss-6-12-IRK');
fprintf('\n\nGauss-6-12-IRK dense output coefficients:\n')
disp(B)

diary off
% end of script