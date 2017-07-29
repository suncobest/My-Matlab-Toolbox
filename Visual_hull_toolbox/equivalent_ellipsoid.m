function [center, semiAxes, vec, inertia, covariance, lamda] = equivalent_ellipsoid(XX)
% Compute the equivalent ellipsoid for points of m-dimension.
%   XX: points of dimension m*npts; (m>=1), assumed as 1D points if XX is vector;
%   center: center of points cloud;
%   semiAxes: semi axes of equivalent ellipsoid in descending order of length;
%   vec: orientation of principal axes correspoinding to semi axes (unit vectors);
%   inertia: moment of inertia with respect to center (matrix);
%   convariance: the convariance matrix wrt center;
%   lamda: the singular value (or eigen value) of convariance matrix

[m,npts] = size(XX);
if npts==1,
    npts = m;
    m = 1;
end;
if npts ==1,
    center = XX;
    semiAxes = 0;
    vec = [];
    inertia = 0;
    covariance = 0;
    lamda = 0;
    return;
end;
center = sum(XX,2)/npts;
YY = XX-center(:,ones(1,npts));
YY = YY*YY';
covariance = YY/npts;
inertia = sum(diag(covariance))*eye(m)-covariance;
[~,S,vec] = svd(YY);
lamda = diag(S)/npts;
semiAxes = sqrt((m+2)*lamda);
    
return;


%% Test
% 对于椭球体x^2/a^2+y^2/b^2+z^2/c^2<=1来说，其体积为4*pi*a*b*c/3；
% x^2对体积积分得到4*pi*a^3*b*c/15；所以5*S(1,1)/Np=a^2。

% 对于椭圆x^2/a^2+y^2/b^2<=1来说，其面积为pi*a*b；
% x^2对面积积分得到Iy=pi*a^3*b/4；所以4*Iy/Np=a^2，Iy为椭圆对x轴的惯性矩

figure(1);hold on;
m = randi([2,3]);
n = 5000;
f = 1.5;
mk = 30;
for i = 1:3,
    render = randn>0;
    show_ellipsoid =  randn>0;
    X_kk=randn(m,n);
    ax0 = sort(randi(5,1,m),'descend');            % ax is also descending sorted
    X_kk=diag(ax0)*X_kk;
    if m==3,
        om = randn(3,1);
        t = randn(3,1)*20;
        X_kk = rodrigues(om)*X_kk+t(:,ones(1,n));
        [ct, ax, V, I] = equivalent_ellipsoid(X_kk);
        
        plot3(X_kk(1,:),X_kk(2,:),X_kk(3,:),'c.')
        plot3(ct(1),ct(2),ct(3),'m*');
        
        % end points
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
            surf(X,Y,Z,'FaceColor', 'b','EdgeColor','g','FaceAlpha',0.3);
            
            if render,
                set(gcf, 'renderer', 'zbuffer','color',[1,1,1]*0.7);
                cameratoolbar('ResetCameraAndSceneLight');
                cameratoolbar('Togglescenelight');
            end;
        end;
    elseif m==2,
        theta = rand*pi;
        R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
        t = randn(2,1)*20;
        X_kk = R*X_kk+t(:,ones(1,n));
        [ct, ax, V, I] = equivalent_ellipsoid(X_kk);
        
        plot(X_kk(1,:),X_kk(2,:),'c.')
        plot(ct(1),ct(2),'m*');
        
        % end points
        cct = ct(:,ones(1,m));
        xy = cct+V*diag(ax)*f;
        plot([cct(1,:);xy(1,:)], [cct(2,:);xy(2,:)],'linewidth',2);
        theta = linspace(0,2*pi,mk);
        xy = ax(:,ones(1,mk)).*[cos(theta); sin(theta)];
        xy = ct(:,ones(1,mk))+V*xy;
        plot(xy(1,:),xy(2,:),'b-');
        axis equal tight off;
    end;
end;
hold off;
