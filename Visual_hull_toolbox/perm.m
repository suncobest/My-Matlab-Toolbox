function  P = perm( x, k )
%   PERM(X,K)从X(:)中任意选取K个不可重复(unrepeatable)元素进行排列：Pnk=Cnk*Pkk
%   不同于对所有元素进行全排列的PERMS。
%   N=numel(x)
%
%   PERM(X,K) creates a matrix with N!/K! rows and K columns containing all 
%   possible permutations of the N elements. 
%
%   Class support for input X:
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%      logical, char
%
%   See also PERMS, NCHOOSEK, RANDPERM, PERMUTE.

if nargin == 1
    k = numel(x);
end

[~,maxsize] = computer;
n = numel(x);

% % Error if n < k
% % if n < k  % 不需要用警告，如果元素不够，那么k个框就空着。
% %     error('There must be enough elements to select! ')  
% % end

% Error if output dimensions are too large
if n>=k && n*factorial(n)/factorial(n-k) > maxsize
    error(message('MATLAB:pmaxsize'))
end

v = x(:); % Make sure v is a column vector

index = permnk( n, k );

P = v(index);

end

%%-----------------------------------------------------------
function  P = permnk( n, k )
% subfunction to help with recursion

if k == 1
    P = (1:n)'; 
    return; 
end

q = permnk( n , k-1);  % recursive calls
m = size(q,1);   % q是m行k-1列
P = zeros(m*(n-k+1), k);  % 初始化，对于n-k+1，当n=1，k=2时等于0，即生成k个框

for i = 1:m
    qi=q(i,:);  % 取出q的第i行
    P((i-1)*(n-k+1)+1:i*(n-k+1), 1:k-1) = qi(ones(n-k+1,1),:);  
    % 将q的每一行复制n-k+1遍再赋给P，其中qi(ones(n-k+1,1),:)等价于repmat(qi,n-k+1,1),但前者计算速度更快。
    
    t=1:n;
    t(q(i,:))=[];  % 从1:n里将q的每行中已经出现过的数字删除，剩下的排在最后
    P((i-1)*(n-k+1)+1:i*(n-k+1), k) = t(:);
end

end
