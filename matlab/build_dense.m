function B = build_dense(method_name)
% Usage: B = build_dense(method_name)
% 
% Returns the dense output table associated with the method given
% by the input argument 'method_name'.  The output table has s+2
% rows (where the method has s stages):
%     B = [a(-1); a(0); a(1); ...; a(s)], 
% where each row a(k) corresponds to the coefficients for the
% basis function 
%     pk(t) = a(k) * Pi(t)
% with 
%     Pi(t) = [1, t, t^2, ..., t^(s+1)]^T
%
% Method types are as specified in butcher.m
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% get the butcher table
Butcher = butcher(method_name);
[m,n] = size(Butcher);

% set the desired precision in digits
dig = 20;

% get the number of stages and extract the 'c' coefficients
s = n-1;
c = Butcher(1:s,1);

% build matrix of interpolation conditions
M = vpa(zeros(s+2,s+2),dig);
M(1,1) = vpa(1,dig);
M(2,:) = vpa(ones(1,s+2),dig);
for i=3:s+2
   for j=2:s+2
      M(i,j) = vpa((j-1)*c(i-2)^(j-2),dig);
   end
end

% solve for coefficients, and place inside return matrix
B = vpa(zeros(s+2,s+2),dig);
for i=1:s+2
   r = vpa(zeros(s+2,1),dig);
   r(i) = vpa(1,dig);
   a = vpa(M\r,dig);
   B(i,:) = vpa(a',20);
end

% end of function
