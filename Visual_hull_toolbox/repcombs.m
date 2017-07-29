function  R=repcombs(x,k)
%% ������
%
% repcomb(x,k)��x(:)������ѡȡk�����ظ�(repeatable)Ԫ�أ�n=numel(x);
% �ֱ�Ž�k�����������n^k�����Σ������������Ρ�
% ÿ��������n��״̬,����x�п������ַ��������󡢻���cell.
%
%   See also REPCOMBNK, COMBNK, PERMS, NCHOOSEK, RANDPERM.


% ���x��cell array�������ó˷�������x(i)*ones(m,1)����Ҫʹ�ø��ƾ�������repmat
% repmat(x,M,N,...)������x�ֱ���ά��1,2,...�ϸ���M,N,...����

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
idx=(1:n)';  % ����x(:)������������ʼ״̬state1��ָ�꣬�Ա�֤repcomb(x,1)=x(:)������(1:n)'Ϊ��������

for j = 1:k-1,  % ��state1��state(k)��Ҫ����k-1�Σ�ע����M>N����M:NΪ�ռ������Ե�k=1ʱ��j=[],�򲻽���ѭ���塣
    m = size(idx,1);  % ���state(j)����m����״̬��
    temp = zeros(n*m,1);    % ��ʼ���м����temp
    for i=1:n,
        temp((i-1)*m+1:i*m) = i*ones(m,1);
    end;
% ��state(j)��state(j+1)
%     idx = [repmat(idx,n,1), temp];  %  �������ҵ������֣��̶��ұߣ��ȵ���ߣ�
    idx = [temp, repmat(idx,n,1)];   %  ��������������֣��̶���ߣ��ȵ��ұߣ�
end;
R = v( idx );
return;
