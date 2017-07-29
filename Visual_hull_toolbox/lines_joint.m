function [X, lamda, foot, stdX] = lines_joint(pts, direc)
% LINE_JOINT compute the joint of lines in the least-squares sense.
%  Input:
%  pts: (ndim*Np), every column is a point on a line. ndim denotes
%  dimension of points. Np is number of lines.
%  direc: (ndim*Np), direction of Np lines.
%  The i-th line: pts(:,i)+lamda(i)*direc(:,i);
%  Output:
%  X: the joint of all lines, which is of dimension ndim*1.
%  lamda and foot: foot(:,i) = pts(:,i)+lamda(:,i)*direc(:,i) is the
%  orthogonal projection of X to the i-th line.
%  stdX is the standard deviation of X.
%
% algorithm: Jean_Yves Bouguet, Visual methods for three-dimensional
% modeling, p58-62 (4.1.2)
%  See also compute_structure.

[ndim, Np] = size(pts);
[m, n] = size(direc);
assert(ndim>1 && Np>1 && m==ndim && n==Np, 'Unexpected inputs!');
nd1 = ones(1,Np);
A = diag(sum(direc.^2))-direc'*direc/Np;
b = sum(pts,2)/Np;
b = sum((b(:,nd1)-pts).*direc);
lamda = b/A;
foot = pts+lamda(ones(ndim,1),:).*direc;
X = sum(foot,2)/Np;
stdX = X(:,nd1)-foot;
stdX = sqrt(sum(sum(stdX.^2))/(Np-1));

return;


%% Test
m = 3;
n = 10;
x = randn(m,1)*10;
dr = randn(m,n);
la = sqrt(sum(dr.^2));
dr = dr./la(ones(m,1),:);
la = randn(1,n)*2;
pt = la(ones(m,1),:).*dr+x(:,ones(1,n))+randn(m,n)/10;
[xx,ll,ff,st] = lines_joint(pt, dr);
err=x-xx
err=la+ll
ptt = ll(ones(m,1),:).*dr*1.2+pt;
if m==2,
    plot(x(1,:),x(2,:),'+',xx(1,:),xx(2,:),'o',[pt(1,:);ptt(1,:)],[pt(2,:);ptt(2,:)],'-',pt(1,:),pt(2,:),'.',ff(1,:),ff(2,:),'*');
elseif m==3,
    plot3(x(1,:),x(2,:),x(3,:),'+',xx(1,:),xx(2,:),xx(3,:),'o',[pt(1,:);ptt(1,:)],[pt(2,:);ptt(2,:)],[pt(3,:);ptt(3,:)],'-',pt(1,:),pt(2,:),pt(3,:),'.',ff(1,:),ff(2,:),ff(3,:),'*');
end;
axis equal;

