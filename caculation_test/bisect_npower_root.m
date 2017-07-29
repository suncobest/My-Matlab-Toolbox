function [Y] = bisect_npower_root(X,N,E)
%BISECT_NPOWER_ROOT get the N-th power root of input argument X.
%   zpf, 2014-01-11
% 
% The function use bisecting method to caculate the N-th power root of X.
% Input argument X and output argument Y can be arrays of Real numbers. The
% power index N must be a Real number. E is the precision of the function.
% If the value of 'X^N'
[i,j]=find(X>=1);
[i1,j1]=find(X<1);
if(N>=1)
   low=1;
   high=X;
   r=(low+high)/2;
   while
    



end

