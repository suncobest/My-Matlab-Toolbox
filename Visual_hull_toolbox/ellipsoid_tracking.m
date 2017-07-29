function [cc_out, ax_out, vec_out, ind] = ellipsoid_tracking(XX, org, ax_in, vec_in, MaxIter, frac, lambda)
% ELLIPSOID_TRACKING optimize orentation of a given ellipsoid to best fit data XX: m-dimensional points cloud.
% Input:
%   XX: points of dimension m*npts; (m>=1), assumed as 1D points if XX is vector;
%   ax_in: initial semi axes of ellipsoid.
%   vec_in:  initial orientation of ellipsoid axes.
%   MaxIter: maximum number of iteration.
%   frac: fraction of long axis to estimate axes direction.
%   org: rotation center of ellipsoid (one end).
%
% Output:
%   cc_out: center of points cloud in the optimized ellipsoid;
%   ax_out: semi axes of points in the segmenting ellipsoid.
%   vec_out: orientation of principal axes correspoinding to semi axes (unit vectors);
%   ind: column logical indices of points in the optimized ellipsoid.
% Note: ellipsoid_tracking is suitable to track rigid body;
% principal_ellipsoid is suitable to track body of small deformation;
% optimal_ellipsoid is suitable to track body of considerable deformation.
%  See also principal_ellipsoid, optimal_ellipsoid.

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

if nargin<7 || isempty(lambda),
    lambda = 1;
end;

if nargin<6 || isempty(frac),
    frac = 0.1;
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

vi = vec_out(:,1);
ctt = cc(:,[1,1])+vi*[1,-1]*ax_out(1)*frac;
YY = YY(:,vi'*(XXo-ctt(:,1)*ones(1,np))<=0 & vi'*(XXo-ctt(:,2)*ones(1,np))>=0);
for ii = 2:m,
    YY = YY-vi*(vi'*YY);
    vc = vec_in(:,ii)-vi*(vi'*vec_in(:,ii));
    Y = vc'*YY;
    [~,I] = min(Y);
    [~,J] = max(Y);
    vi = YY(:,J)-YY(:,I);
    temp = norm(vi);
    if temp~=0,
        vi= vi/temp;
        vec_out(:,ii) = vi;
    else
        % Schmidt orthogonalization
        vc = vec_out(:,ii:m);
        for jj=1:ii-1,
            vc = vc-vec_out(:,jj)*(vec_out(:,jj)'*vc);
        end;
        vc(:,1) = vc(:,1)/norm(vc(:,1));
        vec_out(:,ii:m) = vc;
        for jj =ii+1:m,
            vi = vec_out(:,jj-1);
            vec_out(:,jj:m) = vec_out(:,jj:m)-vi*(vi'*vec_out(:,jj:m));
            vec_out(:,jj) = vec_out(:,jj)/norm(vec_out(:,jj));
        end;
        break;
    end;
end;

for iter = 1:10,
%     fprintf(1,'%d\n',iter);
    % transform to unit sphere and compute points in ellipsoid
    YY = 1./ax_in(:,m1).*vec_out'*(XX-cc_out(:,n1));
    ind = sum(YY.^2,1)<=1;
    np = sum(ind);
    XXo = XX(:,ind);
    cc = sum(XXo,2)/np;
    YY = XXo-cc(:,ones(1,np));
    [~,S,~] = svd(YY*YY');
    ax = sqrt((m+2)*diag(S)/np);
    % compute difference between current and last value
    ex = [cc-cc_out; ax-ax_out];
    cc_out = cc;
    ax_out = ax;
    if norm(ex)<bigeps,
        break;
    end;
end;

vi = vec_out(:,1);
ctt = cc(:,[1,1])+vi*[1,-1]*ax(1)*frac;
YY = YY(:,vi'*(XXo-ctt(:,1)*ones(1,np))<=0 & vi'*(XXo-ctt(:,2)*ones(1,np))>=0);
for ii = 2:m,
    YY = YY-vi*(vi'*YY);
    vc = vec_in(:,ii)-vi*(vi'*vec_in(:,ii));
    Y = vc'*YY;
    [~,I] = min(Y);
    [~,J] = max(Y);
    vi = YY(:,J)-YY(:,I);
    temp = norm(vi);
    if temp~=0,
        vi= vi/temp;
        vec_out(:,ii) = vi;
    else
        % Schmidt orthogonalization
        vc = vec_out(:,ii:m);
        for jj=1:ii-1,
            vc = vc-vec_out(:,jj)*(vec_out(:,jj)'*vc);
        end;
        vc(:,1) = vc(:,1)/norm(vc(:,1));
        vec_out(:,ii:m) = vc;
        for jj =ii+1:m,
            vi = vec_out(:,jj-1);
            vec_out(:,jj:m) = vec_out(:,jj:m)-vi*(vi'*vec_out(:,jj:m));
            vec_out(:,jj) = vec_out(:,jj)/norm(vec_out(:,jj));
        end;
        break;
    end;
end;
    
if flag,
    vec_out = vec_out(:,idx);
    ax_out = ax_out(idx);
end;
return;


Nb = 20;
[Xe, Ye, Ze] = sphere(Nb);
Nb = (Nb+1)^2;
XYZs = [Xe(:)'; Ye(:)'; Ze(:)'];
figure(3);
plot3(XX(1,:)',XX(3,:)',-XX(2,:)', 'g.'); hold on;
plot3(org(1),org(3),-org(2),'+','markersize',30,'linewidth',2);
XYZe = vec*diag(ax_in)*XYZs+cc(:,ones(1,Nb));
Xe(:) = XYZe(1,:);
Ye(:) = XYZe(2,:);
Ze(:) = XYZe(3,:);
surf(Xe,Ze,-Ye,'FaceColor','r', 'EdgeColor','r', 'FaceAlpha',0.05);
axis image off;
hold off;