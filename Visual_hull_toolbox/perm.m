function  P = perm( x, k )
%   PERM(X,K)��X(:)������ѡȡK�������ظ�(unrepeatable)Ԫ�ؽ������У�Pnk=Cnk*Pkk
%   ��ͬ�ڶ�����Ԫ�ؽ���ȫ���е�PERMS��
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
% % if n < k  % ����Ҫ�þ��棬���Ԫ�ز�������ôk����Ϳ��š�
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
m = size(q,1);   % q��m��k-1��
P = zeros(m*(n-k+1), k);  % ��ʼ��������n-k+1����n=1��k=2ʱ����0��������k����

for i = 1:m
    qi=q(i,:);  % ȡ��q�ĵ�i��
    P((i-1)*(n-k+1)+1:i*(n-k+1), 1:k-1) = qi(ones(n-k+1,1),:);  
    % ��q��ÿһ�и���n-k+1���ٸ���P������qi(ones(n-k+1,1),:)�ȼ���repmat(qi,n-k+1,1),��ǰ�߼����ٶȸ��졣
    
    t=1:n;
    t(q(i,:))=[];  % ��1:n�ｫq��ÿ�����Ѿ����ֹ�������ɾ����ʣ�µ��������
    P((i-1)*(n-k+1)+1:i*(n-k+1), k) = t(:);
end

end
