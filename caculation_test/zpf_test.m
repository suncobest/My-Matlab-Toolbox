%% dYdom jacobian:  Y=RX
X = sym('x%d', [3,1],'real');
a = sym('a%d',[3,1],'real');
R = sym('r%d%d',[3,3],'real');
la = norm(a);
Y = R*X;
dYdR = jacobian(Y,R(:));
R = [cos(la)+(1-cos(la))*(a(1)/la)^2, (1-cos(la))*a(1)*a(2)/la.^2-sin(la)*a(3)/la, (1-cos(la))*a(1)*a(3)/la.^2+sin(la)*a(2)/la;
     (1-cos(la))*a(1)*a(2)/la.^2+sin(la)*a(3)/la, cos(la)+(1-cos(la))*(a(2)/la)^2, (1-cos(la))*a(2)*a(3)/la.^2-sin(la)*a(1)/la;
     (1-cos(la))*a(1)*a(3)/la.^2-sin(la)*a(2)/la, (1-cos(la))*a(2)*a(3)/la.^2+sin(la)*a(1)/la, cos(la)+(1-cos(la))*(a(3)/la)^2];
dRda = simplify(jacobian(R(:),a));
dYda = simplify(dYdR*dRda)


%%  track points
tol = 5;  % pixel tolerance of frame difference
filename = 'data.txt';
in = fopen(filename, 'rt');
tline = fgetl(in);
if tline == -1,
    disp('End of file!');
    return;
end;
out = fopen(['new_' filename],'wt');

oldloc = reshape(str2num(tline),2,[]);
m = size(oldloc,2);
fprintf(out,'%-8.3f\t',oldloc(:));
fprintf(out,'\n');

oldid = true(1,m);
% count = 1;
step = zeros(2,m);
while ~feof(in),
    % count = count+1;
    tline = fgetl(in);
    newloc = reshape(str2num(tline),2,[]);
    n = size(newloc,2);
    loc = oldloc+step;
    newid = false(1,n);
    ind = 1:n;   % index number
    dist = sum(abs(permute(repmat(oldloc,[1 1 n]),[3 2 1])-permute(repmat(newloc',[1 1 m]),[1 3 2])), 3);
    ok = dist<tol;
    if any(ok(:)),
        for i = 1:size(newloc,2),
            xi = newloc(:,i);
            [~,id] = min(dist(ok), [], 2);
            loc(:,ind(id)) = xi;
            newid(ind(id)) = true;
            ind(id) = [];
        end;
    end;
    step = zeros(2,n);
    ind = oldid & newid;
    step(:,ind) = loc(:,ind)-oldloc(:,ind);
    for i = 1:n,
        if newid(i),
            fprintf(out,'%-8.3f\t',loc(:,i));
        else
            fprintf(out,'%-8.3f\t%-8.3f\t', -1,-1);
        end;
    end;
    fprintf(out,'\n');
    oldloc = loc;
    oldid = newid;
end
fclose(in);
fclose(out);



%% test of skew matrix
syms a1 a2 a3 b1 b2 b3 c1 c2 c3 real
skew_mat=@(v) [0, -v(3), v(2);v(3), 0, -v(1); -v(2), v(1), 0];
a = [a1;a2;a3];
b = [b1;b2;b3];
c = [c1;c2;c3];
simplify(cross(a,cross(b,c))+cross(b,cross(c,a))+cross(c,cross(a,b)))
simplify(a*b'*skew_mat(c)+skew_mat(c)*b*a'+skew_mat(cross(a,cross(b,c))))

%% Test average of two quaternion
syms q1q2 w1 w2 positive
delta = (w1-w2)^2+4*w1*w2*q1q2^2;
c1sq = 2*(w1*q1q2)^2/(delta+sqrt(delta)*(w2-w1+2*w1*q1q2^2));
c1sq1 = w1*(w1-w2+sqrt(delta))/(sqrt(delta)*(w1+w2+sqrt(delta)));
c2sq = (w2-w1+sqrt(delta))^2/(2*delta+2*sqrt(delta)*(w2-w1+2*w1*q1q2^2));
c2sq1 = w2*(w2-w1+sqrt(delta))/(sqrt(delta)*(w1+w2+sqrt(delta)));
simplify(c1sq-c1sq1)
simplify(c2sq-c2sq1)
simplify(c1sq1+c2sq1+2*sqrt(c1sq1*c2sq1)*q1q2-1)

%% Test [om, domdom1, domdom2]
om1 = randn(3,1);
om2 = randn(3,1);
hand = sign(randn)
dom1 = randn(3,1)/1000;
dom2 = randn(3,1)/1000;
om2t = om2;
if hand~=1,
    om2t(1:2)=-om2t(1:2);
end;
[R1,dR1dom1] = rodrigues(om1);
[R2,dR2dom2] = rodrigues(om2t);
R = R1*R2;
[dRdR1,dRdR2] = dAB(R1,R2);
[om,domdR] = rodrigues(R);
domdom1 = domdR*dRdR1*dR1dom1;
domdom2 = domdR*dRdR2*dR2dom2;
if hand~=1,
    domdom2(:,1:2)=-domdom2(:,1:2);
end;
om2n = om2+dom2;
om2nt = om2n;
if hand~=1,
    om2nt(1:2)=-om2nt(1:2);
end;
omn = rodrigues(rodrigues(om1+dom1)*rodrigues(om2nt))
omm = om + domdom1*dom1+domdom2*dom2
gain = norm(omn-om)/norm(omn-omm)



%% check largest integral point in circle
n=4;
center=zeros(2,n);
radius=zeros(1,n);
maxpts=center;
npt=radius;
upbound=radius;
lowbound=radius;
leftbound=radius;
rightbound=radius;
for i=1:n,
    name=['xdata' num2str(i)];
    eval(['load ' name '.txt']);
    cxy = eval([name '(1,:)']);
    center(:,i) =cxy;
    Mxy = eval([name '(end,:)']);
    maxpts(:,i)=Mxy;
    r = eval([name '(2,1)']);
    radius(i) = r;
    eval([name '=' name '(3:end-1,:);']);
    x = eval([name '(:,1)']);
    y = eval([name '(:,2)']);
    np = length(x);
    npt(i) = np;
    up = ceil(cxy(2)+r);
    upbound(i) = up;
    lo = floor(cxy(2)-r);
    lowbound(i) = lo;
    le = floor(cxy(1)-r);
    leftbound(i) = le;
    ri = ceil(cxy(1)+r);
    rightbound(i) = ri;
    % draw figure(i);
    draw_arc(i,cxy,r,0,360,200);    % draw_arc(AX,center,radius,dstart,dend,ndiv), AX is the figure handle
    hold on;
    for lx=le:ri,
        plot([lx lx],[lo up],'k-');
    end;
    for ly=lo:up,
            plot([le ri],[ly ly],'k-');
    end;
    plot(x,y,'g-',x,y,'co',Mxy(1),Mxy(2),'mo','LineWidth',2.0,'MarkerSize',5.0);

    set(i,'color','w');
    axis equal;
    hold off;
    pause(1);
end;



%% check the speed of cell, structure, and eval defining varibles

clear;
tic;
n = 7000;
a = cell(1,n);
for i=1:n,
    a{i}=randn(3,i);
end;
t1= toc

clear;
tic;
n = 7000;
for i=1:n,
    a(i).data=randn(3,i);
end;
t2= toc

%%
clear;
tic;
n = 7000;
for i=1:n,
    eval(['a' num2str(i) '=randn(3,i);']);
end;
t3= toc



%% squad
syms p1 p2 Q1 Q2 alpha beta gama h1 t real
slerp = @(p,q,a,t) p*sin((1-t)*a)/sin(a)+q*sin(t*a)/sin(a);
SB=slerp(slerp(p1, p2, alpha, t), slerp(Q1, Q2, beta, t), gama, 2*t*(1-t));
SB=simplify(SB)
dxSB=simplify(jacobian(SB,t)/h1);
% dxSB=simplify(gradient(SB,t)/h1);
dxSB0=simplify(subs(dxSB,t,0))
dxSB1=simplify(subs(dxSB,t,1))

%% bezier
syms p1 p2 Q1 Q2 h1 t real
lin = @(p,q,t) (1-t)*p+t*q;
B=lin(lin(p1,p2,t),lin(Q1,Q2,t),2*t*(1-t));
B=simplify(B)
simplify(subs(B,Q1, p1));
simplify(subs(ans,Q2,p2))
dxB=simplify(gradient(B,t)/h1);
dxB0=simplify(subs(dxB,t,0))
dxB1=simplify(subs(dxB,t,1))

%%
syms p1 A1 B1 p2 A2 h1 h2 t real
lin = @(p,q,t) (1-t)*p+t*q;
Bezier=lin(lin(lin(p1,A1,t),lin(A1,B1,t),t),lin(lin(A1,B1,t),lin(B1,p2,t),t),t);
Bezier=simplify(Bezier)
err1 = simplify(Bezier-((1-t)^3*p1+3*(1-t)^2*t*A1+3*(1-t)*t^2*B1+t^3*p2))
err2 = simplify(subs(Bezier,A1,p1)-((1-t)^2*p1+2*(1-t)^2*t*p1+3*(1-t)*t^2*B1+t^3*p2))
err2 = simplify(subs(Bezier,B1,p2)-((3*A1-2*p2-p1)*t^3+3*(p1+p2-2*A1)*t^2+3*(A1-p1)*t+p1))
Bezier1 =[p1, A1, p2+(p2-A2)*h1/h2, p2]*[1 -3 3 -1;0 3 -6 3;0 0 3 -3;0 0 0 1]*[1; t; t^2; t^3];
Bezier1=simplify(Bezier1)
err3 = simplify(Bezier1-((1-t)^3*p1+3*(1-t)^2*t*A1+3*(1-t)*t^2*(p2+(p2-A2)*h1/h2)+t^3*p2))
simplify(subs(Bezier1,A1,p1))
simplify(subs(Bezier1,A2,p2))
simplify(subs(ans,A1,p1))

%%
so = [0.5,1,2,3,5,7,10,50];     % unit: m
f = [50;85;105]/1000;     % unit: m
si = (f*so)./([so; so; so] - f(:,ones(1,length(so))))*1000;       % unit: mm
figure('color','w')
plot(so,si(1,:),'r.-',so,si(2,:),'g.-',so,si(3,:),'b.-')
hold on,axis([1,50,0,120])

%%
I=imread('left01.jpg');

fc_left = Cam_vec(1).fc;
cc_left = Cam_vec(1).cc;
kc_left = Cam_vec(1).kc;
alpha_c_left = Cam_vec(1).alpha_c;

lx=ginput(3)
XL_list = zeros(3,n_cam-1);
xl = lx(1,:)'-1;
for kk = 2:n_cam
    xR = lx(kk,:)'-1;
    om = Cam_vec(kk).om;
    T = Cam_vec(kk).T;
    hand = Cam_vec(kk).hand;
    fc_right = Cam_vec(kk).fc;
    cc_right = Cam_vec(kk).cc;
    kc_right = Cam_vec(kk).kc;
    alpha_c_right = Cam_vec(kk).alpha_c;
    [XL,XR] = stereo_triangulation2(xl,xR,om,T,hand,fc_left,cc_left,kc_left,alpha_c_left,fc_right,cc_right,kc_right,alpha_c_right);
    XL_list(:,kk-1) = XL;   % Xc被重复计算了n_cam-1遍，所以XL_list有n_cam-1列，3*N行
    geom(kk).P0 = XR;
end
geom(1).P0 =  mean(XL_list,2);
XL_list


%% 统一写文字和数字
ex1 = {'a' 1 12 123; 'ab' 4 5 6; 'abc' 7 8 9}
ex2 = cellfun(@ex_func,ex1,'UniformOutput',0);
size_ex2 = cellfun(@length,ex2,'UniformOutput',0);
str_length = max(max(cell2mat(size_ex2)));
ex3 = cellfun(@(x) ex_func2(x,str_length),ex2,'uniformoutput',0)

ex4 = cell2mat(ex3);
fid = fopen('mydata.txt','wt');
for i = 1:size(ex4,1)
    fprintf(fid,'%s   %s   %s   %s\n',ex4(i,1:3),ex4(i,4:6),ex4(i,7:9),ex4(i,10:12));
end
fclose(fid);



%%
a=logical(zeros(200,200,200));
[i,j,k]=ind2sub(size(a),1:numel(a));
i=uint8(i);
j=uint8(j);
k=uint8(k);
a(50:150,50:150,50:150)=1;
sum(sum(sum((a==1))))
tic;b=bwperim(a);toc
tic;b=bwperim(a);toc
tic;b1=bwinperim(a);toc
spy(b(:,:,151))


%%
% a = 10 * randn(3,1000000);
tic

% b = convex_hull(a);
% t1 = toc

c = convex_hull(a);
t2 = toc

%%
% x=1;
% x(2)=1;
% for i=1:25
%     x(end+1)=x(end)+x(end-1);
% end
% x

%%

% syms a b c u v
%
% if(1)
%     error('stop!')
% end
% A=[a,c,u;0,b,v;0,0,1]
% B=A
%
%
% rx=[1,0,0;0,cos(a),-sin(a);0,sin(a),cos(a);]
% ry=[cos(b),0,sin(b);0,1,0;-sin(b),0,cos(b);]
% rz=[cos(c),-sin(c),0;sin(c),cos(c),0;0,0,1;]
%
% % rx=[1,0,0;0,cos(-a),-sin(-a);0,sin(-a),cos(-a);]
% % ry=[cos(-b),0,sin(-b);0,1,0;-sin(-b),0,cos(-b);]
% % rz=[cos(-c),-sin(-c),0;sin(-c),cos(-c),0;0,0,1;]
%
% xyz=rx*ry*rz
% xyz.'
% T=ry*rz

%% 椭球面
m = 30;
[X,Y,Z]=ellipsoid(0,0,0,3,2,1,m);
om = randn(3,1);
t = randn(3,1)*10;
xyz = rodrigues(om)*[X(:)';Y(:)';Z(:)']+t(:,ones(1,(m+1)^2));
X1=X; Y1=Y; Z1=Z;
X1(:) = xyz(1,:);
Y1(:) = xyz(2,:);
Z1(:) = xyz(3,:);
figure(1); hold on;
surf(X,Y,Z,'FaceColor', 'b','EdgeColor','r','FaceAlpha',0.1);
surf(X1,Y1,Z1,'FaceColor', 'b','EdgeColor','r','FaceAlpha',0.1);
axis equal tight vis3d off;
set(gcf, 'renderer', 'zbuffer','color',[1,1,1]*0.7);
cameratoolbar('ResetCameraAndSceneLight');
cameratoolbar('Togglescenelight');


%% 椭球体

X_kk=randn(3,5000);
X_kk=diag([3 2 1])*X_kk;

Np = size(X_kk,2);
X_mean = mean(X_kk,2);
Y = X_kk - (X_mean*ones(1,Np));
YY = Y*Y'
[U,S,V] = svd(YY)
figure,hold on,grid on
plot3(X_kk(1,:),X_kk(2,:),X_kk(3,:),'k.')
axis equal
plot3(X_mean(1),X_mean(2),X_mean(3),'r+')

% 对于椭球体x^2/a^2+y^2/b^2+z^2/c^2<=1来说，其体积为4*pi*a*b*c/3；
% x^2对体积积分得到4*pi*a^3*b*c/15；所以5*S(1,1)/Np=a^2。

% 对于椭圆x^2/a^2+y^2/b^2<=1来说，其面积为pi*a*b；
% x^2对面积积分得到Iy=pi*a^3*b/4；所以4*Iy/Np=a^2，Iy为椭圆对x轴的惯性矩

semi_axes=sqrt(5*diag(S)/Np);

xxx=X_mean+semi_axes(1)*V(:,1);
yyy=X_mean+semi_axes(2)*V(:,2);
zzz=X_mean+semi_axes(3)*V(:,3);

% line([X_mean(1) xxx(1)],[X_mean(2),xxx(2)],[X_mean(3),xxx(3)],'color','r');
% line([X_mean(1) yyy(1)],[X_mean(2),yyy(2)],[X_mean(3),yyy(3)],'color','g');
% line([X_mean(1) zzz(1)],[X_mean(2),zzz(2)],[X_mean(3),zzz(3)],'color','b');

daspect;
arrow3(X_mean',xxx','r',1,2,[]);
arrow3(X_mean',yyy','g',1,2,[]);
arrow3(X_mean',zzz','b',1,2,[]);

% omega=randn(3,1);
% R=rodrigues(omega);
% T=10*rand(3,1);
% X_kk=R*X_kk+T*ones(1,Np);
