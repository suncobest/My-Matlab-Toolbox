function [cc_out, ax_out, vec_out, ind] = ellipsoid_rotation(XX, org, ax_in, vec_in, MaxIter, lambda)
% ELLIPSOID_ROTATION optimize orentation of a given ellipsoid to best fit data XX: m-dimensional points cloud.
% Input:
%   XX: points of dimension m*npts; (m>=1), assumed as 1D points if XX is vector;
%   ax_in: initial semi axes of ellipsoid.
%   vec_in:  initial orientation of ellipsoid axes.
%   MaxIter: maximum number of iteration.
%   org: rotation center of ellipsoid (one end).
%
% Output:
%   cc_out: center of points cloud in the optimized ellipsoid;
%   ax_out: semi axes of points in the segmenting ellipsoid.
%   vec_out: orientation of principal axes correspoinding to semi axes (unit vectors);
%   ind: column logical indices of points in the optimized ellipsoid.
%  See also ellipsoid_rotation2.

bigeps = 1e-8;

[m,npts] = size(XX);
if npts==1,
    npts = m;
    m = 1;
end;
if npts ==1,
    cc_out = XX;
    vec_out = [];
    ind = true;
    return;
end;

if nargin<6 || isempty(lambda),
    lambda = 1;
end;

if nargin<5 || isempty(MaxIter),
    MaxIter = 20;
else
    assert(MaxIter>=1,'The function need at least iterate once!');
end;

if nargin<4 || isempty(vec_in),
    vec_in = eye(m);
else
    assert(ismatrix(vec_in) && isequal(size(vec_in),[m,m]), 'Unexpected dimension of initial axes orientation!');
    [vec_out,~,vec_in] = svd(vec_in);
    vec_in = vec_out*vec_in';
end;

m1 = ones(1,m);
n1 = ones(1,npts);
flag = 0;     % switch to turn on axis refinement
if nargin<3 || isempty(ax_in),
    cc_out = sum(XX,2)/npts;
    YY = XX-cc_out(:,n1);
    [~,S,vec_in] = svd(YY*YY');
    ax_in = sqrt((m+2)*diag(S)/npts);
else
    if isscalar(ax_in),
        ax_in = ax_in*m1';
    else
        flag = 1;
        assert(isvector(ax_in) && length(ax_in)==m, 'Unexpected dimension of initial ellipsoid semi axes!');
        [ax_in,idx] = sort(ax_in(:),'descend');
        vec_in = vec_in(:,idx);
        [~,idx] = sort(idx);
    end;
end;

assert(isvector(org) && length(org)==m, 'Unexpected dimension of initial center position!');
org = org(:);

cc_out = org+vec_in(:,1)*ax_in(1)*lambda;
ax = ax_in;
vec_out = vec_in;
volume = prod(ax_in);
for iter = 1:MaxIter,
    %     fprintf(1,'%d\n',iter);
    % transform to unit sphere and compute points in ellipsoid
    YY = 1./ax(:,m1).*vec_out'*(XX-cc_out(:,n1));
    ind = sum(YY.^2,1)<=1;
    np = sum(ind);
    if np==0,
        cc_out = [];
        ax_out = [];
        vec_out = [];
        ind = [];
        return;
    end;
    XXo = XX(:,ind);
    cc = sum(XXo,2)/np;
    YY = XXo-cc(:,ones(1,np));
    [~,S,vec] = svd(YY*YY');
    ax_out = sqrt((m+2)*diag(S)/np);
    axs = ax_out*(volume/prod(ax_out))^(1/m);
    % optimize orientation using max projection criterion
    [~,J] = max(vec_in(:,1)'*YY);
    vi = XXo(:,J)-org;
    vec(:,1) = vi/norm(vi);
    % Schmidt orthogonalization
    for jj =2:m,
        vi = vec(:,jj-1);
        vec(:,jj:m) = vec(:,jj:m)-vi*(vi'*vec(:,jj:m));
        vec(:,jj) = vec(:,jj)/norm(vec(:,jj));
    end;
    cc = org+vec(:,1)*ax_in(1)*lambda;
    ex = abs(sum(vec_out.*vec,1))-m1;
    ex = [ex(:); cc-cc_out; axs-ax];
    vec_out = vec;
    cc_out = cc;
    ax = axs;
    if norm(ex)<bigeps,
        break;
    end;
end;

id = sum(vec_out.*vec_in)<0;
vec_out(:,id) = -vec_out(:,id);
if flag,
    vec_out = vec_out(:,idx);
    ax_out = ax_out(idx);
end;
return;