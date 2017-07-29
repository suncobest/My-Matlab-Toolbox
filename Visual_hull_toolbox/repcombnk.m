function  R = repcombnk( x, k )
%% 密码轮
%
% 与函数repcomb作用一样，不同的是repcomb使用循环算法，repcombnk使用递归算法。
% repcombnk(x,k)从x(:)中任意选取k个可重复(repeatable)元素，
% 分别放进k个格子里，共有numel(x)^k种情形，给出所有情形。
% 每个格子有n种状态,其中x必须是字符串或数值矩阵，其中n = numel(x)
% 
%   See also REPCOMBS, COMBS, COMBNK, PERMS, NCHOOSEK.

% repmat(x,M,N,...)将数组x分别在维度1,2,...上复制M,N,...个。

if nargin<2
    k=numel(x);
end

[~,maxsize] = computer;
n = numel(x);

% Error if output dimensions are too large
if  n*n^k > maxsize
    error(message('MATLAB:pmaxsize'))
end

v = x(:); % Make sure v is a column vector

idx = repcom( n, k );

R = v(idx);

end

%%------------------------------------------------------------------------
function  R = repcom( n, k )
if k==1
    R = (1:n)';  % 终点 R = repcom( n, 1 ) = (1:n)'
    return;
end

s = repcom( n, k-1 );   % recursive calls
m = size(s,1);

temp = zeros(n*m,1);    % 初始化中间变量temp
for i=1:n
    temp((i-1)*m+1:i*m) = i*ones(m,1);
end

%     R = [repmat(s,n,1),temp];  %  从左向右调密码轮（固定右边，先调左边）

R = [temp,repmat(s,n,1)];   %  从右向左调密码轮（固定左边，先调右边）

end