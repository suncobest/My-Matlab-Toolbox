function out = combs(in,k)
%COMB All combinations of the N elements in V taken K at a time.
%   OUT = COMBS(IN,K) produces a matrix, with K columns. Each row of OUT has
%   K of the elements in the vector IN. OUT has N!/K!(N-K)! rows.  K must be
%   a nonnegative integer.
% See also perms.

if nargin<2,
    k = numel(in);
end;
in = in(:).';       % in is a row vector
n = length(in);
assert(k<=n && k>0, 'Unexpted value of the 2nd varible!');

if n == k,                        %  ��n������ȡn������ΪCn(n)
   out = in;
elseif n == k+1,             %  ��n������ȡn-1������ΪCn(n-1)
   out = in(ones(n,1),:);     % ��in����n�У��ȼ���out = repmat(in,n,1); ǰ������Ͽ졣
   out(1:n+1:n*n) = [];
   % ��out�е�һ��Ԫ�ؿ�ʼ��ÿ��n+1��Ԫ��ɾ��һ������ɾ������a�����Խ���Ԫ��,��ӵ�1
   % �е��� n�зֱ�ɾ����1:n��ʣ��n(n-1)��Ԫ��������������Ӿ�����ɾ��Ԫ�غ󣬾���
   % �ͱ������������
   out = reshape(out,n,n-1);
elseif k == 1,          % Cn(1)
   out = in.';
else                              % this algorithm need large space to store data
   rows = 2.^n;             % n bit will have 2^n status, 1 denote chosen, 0 not
   ncycles = rows;         % list all possible status, add chose those row whose sum is k 
   logid = false(rows,n);
   for cnt = 1:n,
      status = [false,true];        % ~[1;0];
      ncycles = ncycles/2;
      nreps = rows./(2*ncycles);
      status = status(ones(1,nreps),:);
      status = status(:);
      status = status(:,ones(1,ncycles));
      logid(:,n-cnt+1) = status(:);
   end;
   logid = logid(sum(logid,2)==k, :);
   nr = size(logid,1);
   if n<20,
       out = in(ones(nr,1), :)';                    % ��ά�߼�ָ�����ûὫ��ά������������
       out = reshape(out(logid'),k,nr)';
   else
       out = zeros(nr,k);
       for cnt =1:nr,
           out(cnt,:) = in(logid(cnt,:));
       end;
   end;
end;
   
return;

