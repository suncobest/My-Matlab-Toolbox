function  R=repcombs(x,k)
%% 密码轮
%
% repcomb(x,k)从x(:)中任意选取k个可重复(repeatable)元素，n=numel(x);
% 分别放进k个格子里，共有n^k种情形，给出所有情形。
% 每个格子有n种状态,其中x有可能是字符串、矩阵、或者cell.
%
%   See also REPCOMBNK, COMBNK, PERMS, NCHOOSEK, RANDPERM.


% 如果x是cell array，则不能用乘法运算如x(i)*ones(m,1)，需要使用复制矩阵命令repmat
% repmat(x,M,N,...)将数组x分别在维度1,2,...上复制M,N,...个。

if nargin<2,
    k = numel(x);
end;
[~,maxsize] = computer;
n = numel(x);
% Error if output dimensions are too large
if  n*n^k > maxsize,
    error(message('MATLAB:pmaxsize'));
end;

v = x(:); % Make sure v is a column vector
idx=(1:n)';  % 建立x(:)的索引，即初始状态state1的指标，以保证repcomb(x,1)=x(:)，其中(1:n)'为列向量。

for j = 1:k-1,  % 从state1到state(k)需要迭代k-1次；注意若M>N，则M:N为空集。所以当k=1时，j=[],则不进入循环体。
    m = size(idx,1);  % 求出state(j)行数m，即状态数
    temp = zeros(n*m,1);    % 初始化中间变量temp
    for i=1:n,
        temp((i-1)*m+1:i*m) = i*ones(m,1);
    end;
% 从state(j)到state(j+1)
%     idx = [repmat(idx,n,1), temp];  %  从左向右调密码轮（固定右边，先调左边）
    idx = [temp, repmat(idx,n,1)];   %  从右向左调密码轮（固定左边，先调右边）
end;
R = v( idx );
return;
