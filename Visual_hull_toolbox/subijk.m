function subpoints = subijk(IJK,ndabc)
% IJK ----- (3,n);   ndabc ----- scalar or vector length of 3;
% ndabcΪĸ�㵽�Ӳ��ϸ������
% �����±�ΪIJK��ĸ����У�ʹ��ϸ��ndabc�󣬵õ�ĸ�����а����������ӵ���±�subpoints

% ʹ���������ʱ������ȷ��IJK��uint8���͵�����(ת�����׳���)

[m,n] = size(ndabc);
if min(m,n)~=1
    error('Unexpected subdivision levels!')
elseif m==n    % m==n==1
    ndabc = ndabc * ones(3,1);
elseif max(m,n)~=3
     error('Unexpected subdivision levels!')
end
 
[m,n] = size(IJK);  % n is number of points in parent layer; nΪĸ��ĵ���
if m~=3
    error('Points have to be 3 dimensional! (3 rows)');
end

% if max(ndabc)>255
%     warning('Subdivision levels exceed 255; cut by funtion uint8!');
% end

% ndabc = uint8(ndabc); % һ���ܱ�֤ϸ����������255������໮�ֳ�255^3=16581375���㣩��ʵ���ϴﲻ��������ͻ�out of memory

nda = ndabc(1);
ndb = ndabc(2);
ndc = ndabc(3);

N = nda*ndb*ndc;  % N is number of child points in one parent cell;
IJK = reshape(IJK,3,1,[]);  % put points dimention into pages
I = repmat((IJK(1,1,:)-1) * nda, 1, nda) + repmat((1:nda),1,1,n);   % ��ÿ�������ͳһ������ÿһҳ���õ�(I-1)*nda+1 : I*nda
J = repmat((IJK(2,1,:)-1) * ndb, 1, ndb) + repmat((1:ndb),1,1,n);
K = repmat((IJK(3,1,:)-1) * ndc, 1, ndc) + repmat((1:ndc),1,1,n);

i = repmat(I,1,ndb*ndc);
j = reshape(repmat(J,nda,ndc),1,N,[]);
k = reshape(repmat(K,nda*ndb,1),1,N,[]);

subpoints = reshape([i;j;k],3,[]);

return;

%% ����ĸ����ӵ��㷨
% I = IJK(1);
% J = IJK(2);
% K = IJK(3);
% nda = ndabc(1);
% ndb = ndabc(2);
% ndc = ndabc(3);

% i = repmat( (I-1)*nda+1 : I*nda, 1, ndb*ndc);
% j = reshape(repmat( (J-1)*ndb+1 : J*ndb, nda, ndc), 1,[]);
% k = reshape(repmat( (K-1)*ndc+1 : K*ndc, nda*ndb, 1), 1,[]);
% subpoints = [i;j;k];
