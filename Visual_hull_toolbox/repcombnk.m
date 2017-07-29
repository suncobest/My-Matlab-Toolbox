function  R = repcombnk( x, k )
%% ������
%
% �뺯��repcomb����һ������ͬ����repcombʹ��ѭ���㷨��repcombnkʹ�õݹ��㷨��
% repcombnk(x,k)��x(:)������ѡȡk�����ظ�(repeatable)Ԫ�أ�
% �ֱ�Ž�k�����������numel(x)^k�����Σ������������Ρ�
% ÿ��������n��״̬,����x�������ַ�������ֵ��������n = numel(x)
% 
%   See also REPCOMBS, COMBS, COMBNK, PERMS, NCHOOSEK.

% repmat(x,M,N,...)������x�ֱ���ά��1,2,...�ϸ���M,N,...����

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
    R = (1:n)';  % �յ� R = repcom( n, 1 ) = (1:n)'
    return;
end

s = repcom( n, k-1 );   % recursive calls
m = size(s,1);

temp = zeros(n*m,1);    % ��ʼ���м����temp
for i=1:n
    temp((i-1)*m+1:i*m) = i*ones(m,1);
end

%     R = [repmat(s,n,1),temp];  %  �������ҵ������֣��̶��ұߣ��ȵ���ߣ�

R = [temp,repmat(s,n,1)];   %  ��������������֣��̶���ߣ��ȵ��ұߣ�

end