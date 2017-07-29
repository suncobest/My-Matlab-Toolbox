function voxels = subvoxcenter(center,ndabc,abcsize,omlist)

% center,omlist ----- (3,N); or omlist ----- (3,1)
% abcsize,ndabc ------(3,N) or (3,1) or (1,1)
% voxels ------- ( 3, sum(nda*ndb*ndc) )
% N is the total number of center seeds to generate grids; center are the coordinates of N seeds (centers of cuboids); 
% abcsize are the sizes of N cuboids (prarent grids); ndabc are the subdivision levels of N cuboids (grids);
% omlist are the directions of N cuboids(grids); There are total sum(nx*ny*nz) voxel points after subdivision.
%
% Example:
% voxels = subvoxcenter(center), return center itself
% subvoxcenter(center,ndabc)��ʾ�����������ϸ��
% ���ڲ�ͬλ�ã���ͬ��С���Լ���ͬ�����ĸ��Ԫ�����������ô˺�������ϸ�֣��õ�����������

if isempty(center),
    voxels = [];
    return;
end;

[dims,N] = size(center);
if dims ~= 3,
    error('Points have to be 3 dimensional! (3 rows)');
end;

if nargin==1,
    voxels = center;
    return;
end;

[m,n] = size(ndabc);
if m==1,
    ndabc = ndabc(ones(3,1),:);
elseif  m~=3,
    error('Subdivision levels have to be 1 or 3 dimensional!');
end;

sub_center = cell(1,N);
% ����repmat������ռ���ڴ棬��ȻҪ����жϡ�
if nargin == 2,
    if n==1,
        for kk = 1:N,
            sub_center{kk} = subdivide(center(:,kk),ndabc,1,0);
        end;
    elseif n==N,
        for kk = 1:N,
            sub_center{kk} = subdivide(center(:,kk),ndabc(:,kk),1,0);
        end;
    else
        error('Points number unmatched!');
    end;
    voxels = cell2mat(sub_center);
    return;
end;


[m1,n1] = size(abcsize);
if m1~=1 && m1~=3,
    error('Size of cubes have to be 1 or 3 dimensional!');
end;
% ��abcsize��һά�����踴�Ƹ��Ƴ����С����Ӻ���subdivide

if nargin == 3,
    if n1==1,
        if n==1,
            for kk = 1:N,
                sub_center{kk} = subdivide(center(:,kk),ndabc,abcsize,0);
            end;
        elseif n==N,
            for kk = 1:N,
                sub_center{kk} = subdivide(center(:,kk),ndabc(:,kk),abcsize,0);
            end;
        else
            error('Points number unmatched!');
        end;
    elseif n1==N,
        if n==1,
            for kk = 1:N,
                sub_center{kk} = subdivide(center(:,kk),ndabc,abcsize(:,kk),0);
            end;
        elseif n==N,
            for kk = 1:N,
                sub_center{kk} = subdivide(center(:,kk),ndabc(:,kk),abcsize(:,kk),0);
            end
        else
            error('Points number unmatched!');
        end
    else
        error('Points number unmatched!');
    end;
    voxels = cell2mat(sub_center);
    return;
end;
    

[m2,n2] = size(omlist);
if m2 ~= 3,
    error('Omega (rotation) have to be 3 dimensional! (3 rows)');
end;

if n2==1,
    if n1==1,
        if n==1,
            for kk = 1:N,
                sub_center{kk} = subdivide(center(:,kk),ndabc,abcsize,omlist);
            end;
        elseif n==N,
            for kk = 1:N,
                sub_center{kk} = subdivide(center(:,kk),ndabc(:,kk),abcsize,omlist);
            end;
        else
            error('Points number unmatched!');
        end;
    elseif n1==N,
        if n==1
            for kk = 1:N,
                sub_center{kk} = subdivide(center(:,kk),ndabc,abcsize(:,kk),omlist);
            end;
        elseif n==N,
            for kk = 1:N,
                sub_center{kk} = subdivide(center(:,kk),ndabc(:,kk),abcsize(:,kk),omlist);
            end;
        else
            error('Points number unmatched!');
        end;
    else
        error('Points number unmatched!');
    end;
    
elseif n2==N,
    if n1==1,
        if n==1,
            for kk = 1:N,
                sub_center{kk} = subdivide(center(:,kk),ndabc,abcsize,omlist(:,kk));
            end;
        elseif n==N,
            for kk = 1:N,
                sub_center{kk} = subdivide(center(:,kk),ndabc(:,kk),abcsize,omlist(:,kk));
            end;
        else
            error('Points number unmatched!');
        end
    elseif n1==N,
        if n==1,
            for kk = 1:N,
                sub_center{kk} = subdivide(center(:,kk),ndabc,abcsize(:,kk),omlist(:,kk));
            end;
        elseif n==N,
            for kk = 1:N,
                sub_center{kk} = subdivide(center(:,kk),ndabc(:,kk),abcsize(:,kk),omlist(:,kk));
            end;
        else
            error('Points number unmatched!');
        end;
    else
        error('Points number unmatched!');
    end;
    
else
    error('Points number unmatched!');
end;

voxels = cell2mat(sub_center);

end



% ����һ������������nxyzΪx��y��z�����ϸ������abcΪ��������ĳ߶ȣ�omΪ������ķ���CΪ���������ĵ�

function points = subdivide(C,nxyz,abc,om)

nx = nxyz(1);
ny = nxyz(2);
nz = nxyz(3);

% nx = double(nxyz(1));
% ny = double(nxyz(2));
% nz = double(nxyz(3));

Np = nx*ny*nz;  % number of voxels
[i,j,k] = ind2sub([nx,ny,nz], 1:Np);   %����1:(nx*ny*nz)��������������i,j,kҲȫ��������


% ��ʼ������������λ��ԭ�㣬ϸ�ֵ�λ������(initialization)
% ��ָ��i��j��k�������������ص����ĵ�(x,y,z)
x = (i-1/2)/nx -1/2;
y = (j-1/2)/ny -1/2;
z = (k-1/2)/nz -1/2;

% resize grids
if length(abc)==1,
    points = abc * [x;y;z];     % vertcat(x,y,z);
else
    points = diag(abc) * [x;y;z];
end;

% rotate grids
if any(om),
    points = rodrigues(om) * points;
end;

% translate grids
if any(C),
    points = points + C(:,ones(1,Np));
end;
end