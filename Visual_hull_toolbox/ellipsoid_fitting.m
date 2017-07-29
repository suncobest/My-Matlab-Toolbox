function [cc_out, ax_out, vec_out, ind] = ellipsoid_fitting(XX, cc_in, ax_in, vec_in, MaxIter, SW)
% ELLIPSOID_FITTING optimize position and orentation of a given ellipsoid
% to best fit data XX: m-dimensional points cloud.
% Input:
%   XX: points of dimension m*npts; (m>=1), assumed as 1D points if XX is vector;
%   cc_in: initial center position of ellipsoid to be optimized.
%   ax_in: ellipsoid size: semi axes of ellipsoid.
%   vec_in: initial orientation of ellipsoid axes.
%   MaxIter:  maximum number of iteration.
%   SW: switch to turn on/off update of orientation. (default: on)
%
% Output:
%   cc_out: center of points cloud in the optimized ellipsoid;
%   ax_out: semi axes of points in the segmenting ellipsoid.
%   vec_out: orientation of principal axes correspoinding to semi axes (unit vectors);
%   ind: column logical indices of points in the optimized ellipsoid.
% Note: ellipsoid_fitting is suitable to track rigid body;
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
    ax_out = 0;
    vec_out = [];
    ind = true;
    return;
end;

if nargin<6 || isempty(SW),
    SW = 1;
else
    SW = ~~SW;
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

if nargin<2 || isempty(cc_in),
    cc_in = sum(XX,2)/npts;
else
    assert(isvector(cc_in) && length(cc_in)==m, 'Unexpected dimension of initial center position!');
    cc_in = cc_in(:);
end;

cc_out = cc_in;
ax_out = ax_in;
vec_out = vec_in;
for iter = 1:MaxIter,
%     fprintf(1,'%d\n',iter);
    % transform to unit sphere and compute points in ellipsoid
    YY = 1./ax_in(:,m1).*vec_out'*(XX-cc_out(:,n1));
    ind = sum(YY.^2,1)<=1;
    np = sum(ind);
    XXo = XX(:,ind);
    cc = sum(XXo,2)/np;
    YY = XXo-cc(:,ones(1,np));
    [~,S,vec] = svd(YY*YY');
    ax = sqrt((m+2)*diag(S)/np);
    % compute difference between current and last value
    if SW, % optimize orientation using max projection criterion
        if vec(:,1)'*vec_in(:,1)<0,
            vec(:,1) = -vec(:,1);
        end;
        vi = vec(:,1);
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
                vec(:,ii) = vi;
            else
                % Schmidt orthogonalization
                vc = vec(:,ii:m);
                for jj=1:ii-1,
                    vc = vc-vec(:,jj)*(vec(:,jj)'*vc);
                end;
                vc(:,1) = vc(:,1)/norm(vc(:,1));
                vec(:,ii:m) = vc;
                for jj =ii+1:m,
                    vi = vec(:,jj-1);
                    vec(:,jj:m) = vec(:,jj:m)-vi*(vi'*vec(:,jj:m));
                    vec(:,jj) = vec(:,jj)/norm(vec(:,jj));
                end;
                break;
            end;
        end;
        ex = abs(sum(vec_out.*vec,1))-m1;
        ex = [ex(:); cc-cc_out; ax-ax_out];
        vec_out = vec;
    else
        ex = [cc-cc_out; ax-ax_out];
    end;
    cc_out = cc;
    ax_out = ax;
    if norm(ex)<bigeps,
        break;
    end;
end;

if flag,
    vec_out = vec_out(:,idx);
    ax_out = ax_out(idx);
end;
return;



%% Test

m = randi([2,3]);
n = 30;
count = 2;
mk = 30;
show_ellipsoid = randn>0;
t0 = randi([5,8])*(1+rand(m,count));
ax0 = sort(randi(4,m,count),'descend');            % ax is also descending sorted

if m==3,
    f = 1.5;
    [xg,yg,zg]=meshgrid(-1:2/n:1);
    X0 = [xg(:)';yg(:)';zg(:)'];
    id = sum(X0.^2,1)<=1;
    X0 = X0(:,id);
    N = size(X0,2);
    om = randn(3,count);
    X_kk = zeros(m,N*count);
    for i =1:count,
        t = t0(:,i);
        Y_kk= rodrigues(om(:,i))*diag(ax0(:,i))*X0+t(:,ones(1,N));
        X_kk(:,(i-1)*N+1:i*N) = Y_kk;
    end;
    %
    id = round(rand*count*N);
    cc = X_kk(:,id);
    [ct, ax, V, id] = ellipsoid_fitting(X_kk,cc,ax0(:,randi(count))*(1+randn/10),[],100);
    nid = ~id;
    close; figure(1);hold on;
    plot3(X_kk(1,id),X_kk(2,id),X_kk(3,id),'c.');
    plot3(X_kk(1,nid),X_kk(2,nid),X_kk(3,nid),'m.')
%     plot3(ct(1),ct(2),ct(3),'m*');
    
    cct = ct(:,ones(1,3));
    xyz = cct+V*diag(ax)*f;
    plot3([cct(1,:);xyz(1,:)], [cct(2,:);xyz(2,:)],[cct(3,:);xyz(3,:)],'linewidth',2);
    %     arrow3(cct',xyz','x',1,2);
    
    axis equal tight vis3d off;
    if show_ellipsoid,
        [X,Y,Z]=ellipsoid(0,0,0,ax(1),ax(2),ax(3),mk);
        xyz = V*[X(:)';Y(:)';Z(:)']+ct(:,ones(1,(mk+1)^2));
        X(:) = xyz(1,:);
        Y(:) = xyz(2,:);
        Z(:) = xyz(3,:);
        surf(X,Y,Z,'FaceColor', 'b','EdgeColor','g','FaceAlpha',0.1);
    end;
elseif m==2,
    f = 1;
     [xg,yg]=meshgrid(-1:2/n:1);
    X0 = [xg(:)';yg(:)'];
    id = sum(X0.^2,1)<=1;
    X0 = X0(:,id);
    N = size(X0,2);
    theta = sign(randn)*pi/2*rand(1,count);
    X_kk = zeros(m,N*count);
    
    for i = 1:count,
        t = theta(i);
        R = [cos(t), -sin(t); sin(t), cos(t)];
        t = t0(:,i);
        Y_kk=R*diag(ax0(:,i))*X0+t(:,ones(1,N));
        X_kk(:,(i-1)*N+1:i*N) = Y_kk;
    end;
    %
    id = round(rand*count*N);
    cc = X_kk(:,id);
    [ct, ax, V, id] = ellipsoid_fitting(X_kk, cc, ax0(:,randi(count))*(1+randn/10),[],50);
    nid = ~id;
    
    close; figure(1);hold on;
    plot(X_kk(1,id),X_kk(2,id),'c.');
    plot(X_kk(1,nid),X_kk(2,nid),'g.')
    plot(ct(1),ct(2),'m*');
    
    cct = ct(:,ones(1,m));
    xy = cct+V*diag(ax)*f;
    plot([cct(1,:);xy(1,:)], [cct(2,:);xy(2,:)],'linewidth',2);
    % draw ellipse
    t = linspace(0,2*pi,mk);
    xy = ax(:,ones(1,mk)).*[cos(t); sin(t)];
    xy = ct(:,ones(1,mk))+V*xy;
    plot(xy(1,:),xy(2,:),'b-');
    axis equal tight off;
end;

