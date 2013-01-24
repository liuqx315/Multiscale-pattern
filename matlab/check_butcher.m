function [q,p,Bs,As,Ls] = check_butcher(B)
% Usage: [q,p,Bs,As,Ls] = check_butcher(B)
% 
% Checks the Butcher table B to determine the analytical order of
% accuracy for the method (q) and its embedding (p), whether the
% method is B stable, and estimates of whether the method is A
% and/or L stable.  It is assumed that B has block structure
%     B = [c, A; 0, b]
% for a standard Runge-Kutta method, or 
%     B = [c, A; 0, b; 0, b2]
% if the method has an embedded error indicator.  
%
% If the method has no embedding, then we set p=0.
%
% Daniel R. Reynolds
% SMU Mathematics
% May 2012
% All Rights Reserved

% set tolerance on 'equality'
tol = 1e-8;

% extract components of Butcher table
[m,n] = size(B);
if (m == n)         % no embedding
   b2 = 0;
elseif (m == n+1)  % has an embedding
   b2 = B(m,2:n)';
else   % illegal input
   error('illegal Butcher table input')
end
s = n-1;
c = B(1:s,1);
b = B(s+1,2:n)';
A = B(1:s,2:n);

% determine P, Q and R for the Butcher simplifying assumptions
%   B(P):
P = 0;
for i=1:1000
   LHS = b'*(c.^(i-1));
   RHS = 1/i;
   if (abs(RHS-LHS)>tol)
      break;  
   end
   P = P+1;
end

%   C(Q):
Q = 0;
for k=1:1000
   alltrue = 1;
   for i=1:s
      LHS = A(i,:)*(c.^(k-1));
      RHS = c(i)^k/k;
      if (abs(RHS-LHS)>tol)
         alltrue=0;
         break;
      end
   end
   if (alltrue == 1)
      Q = Q+1;
   else
      break;
   end
end

%   D(R):
R = 0;
for k=1:1000
   alltrue = 1;
   for j=1:s
      LHS = 0;
      for i=1:s
         LHS = LHS + A(i,j)*b(i)*c(i)^(k-1);
      end
      RHS = b(j)/k*(1-c(j)^k);
      if (abs(RHS-LHS)>tol)
         alltrue = 0;
         break;  
      end
   end
   if (alltrue == 1)
      R = R+1;
   else
      break;
   end
end

% determine q
q = 0;
for i=1:P
   if ((q > Q+R+1) | (q > 2*Q+2)),  
      break;  
   end
   q = q+1;
end


% if there's an embedding, determine the order
if (length(b2) > 1) 
   P = 0;
   for i=1:1000
      LHS = b2'*(c.^(i-1));
      RHS = 1/i;
      if (abs(RHS-LHS)>tol),  break;  end
      P = P+1;
   end
   R = 0;
   for k=1:1000
      alltrue = 1;
      for j=1:s
         LHS = 0;
         for i=1:s
            LHS = LHS + A(i,j)*b2(i)*c(i)^(k-1);
         end
         RHS = b2(j)/k*(1-c(j)^k);
         if (abs(RHS-LHS)>tol)
            alltrue = 0;
            break;  
         end
      end
      if (alltrue == 1)
         R = R+1;
      else
         break;
      end
   end
   p = 0;
   for i=1:P
      if ((p > Q+R+1) | (p > 2*Q+2)),  break;  end
      p = p+1;
   end
  
else
   p = 0;
end



% determine B stability
M = zeros(s,s);
for j=1:s
   for i=1:s
      M(i,j) = b(i)*A(i,j) + b(j)*A(j,i) - b(i)*b(j);
   end
end
lam = eig(M);
Bs = 1;
for i=1:s
   if (lam(i) < 0)
      Bs = 0;
      break;
   end
end


% estimate A stability
i = sqrt(-1);
ztests = [-eps+1000*i, -1000+i, -1000-1000*i, -10000+2+i, -100000];
As = 1;
e = ones(s,1);
I = eye(s);
for eta = ztests
   Reta = 1 + eta*b'*((I-eta*A)\e);
   if (abs(Reta)>1)
      As = 0;
      break;
   end
end


% estimate L stability
ztests = -logspace(0,8,9);
vals = zeros(size(ztests));
for i = 1:length(ztests); 
   eta = ztests(i);
   Reta = 1 + eta*b'*((I-eta*A)\e);
   vals(i) = abs(Reta);
end
Ls = 1;
for i=2:length(ztests)
   if (vals(i) > vals(i-1))
      Ls = 0;
      break;
   end
end



% end of function
