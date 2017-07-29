function out = bezier_interp(in,ndiv) 
% BEZIER_INTERP -  bezier curve interpolation.
% This routine use De Casteljau algorithm to interplate between the first
% and last points of IN which include all nodal points for bezier curve. 
% Note: only the first and last points are on the bezier curve.
%
% IN    : ndim by Np, the input nodal points for interpolation.
% NDIV  : 1 by 1 interger, the subdivision number of the curve.
% OUT   : ndim by ndiv+1, the output curve points including first and last nodal points.
%
% If IN is a matrix, then each column is a point, and each row is a dimension.
% The order of bezier curve is Np-1  

% by zpf, form BIT, 2015-5-23

if nargin < 2,
   ndiv = 10;
   fprintf(1,'\nUsing default subdivision number for interpolation: 10;\n\n');
end;

if ndiv<0 || floor(ndiv)~=ndiv,
    error('Unexpected subdivision number!');
end;

[m,Np] = size(in);
ndim = m;

if Np==1,
    if m==1,
        out = in;
        return;
    end;
    in=in'; Np=m; ndim=1;
end;


% order = Np -1;

out = zeros(ndim,ndiv+1);
out(:,1)=in(:,1);

dt = 1/ndiv;
t = 0;
for ind = 2:ndiv,
    t = t+dt;
    xk = in;
    for k = Np:-1:2,                                      % De Casteljau algorithm
        xk = (1-t)*xk(:,1:k-1)+t*xk(:,2:k);
    end;
    out(:,ind) = xk; 
end;

out(:,end) = in(:,end);

if ndim==1 && m==Np,
    out = out';
end;

return;



%% Test

N = 6;
Ndim = 3;
a=10*rand(Ndim,N);

min_a = min(a');
max_a = max(a');
delta = max(abs(a(:)))/100;

ndiv = 30;
b=bezier_interp(a,ndiv);

figure(1);
hold off;

if Ndim==2,
    plot(b(1,:),b(2,:),'g-','linewidth',3.0);     
elseif Ndim==3,
    plot3(b(1,:),b(2,:),b(3,:),'g-','linewidth',3.0);    
end;

hold on;

for i=1:N,
    if Ndim==2,
        plot(a(1,i),a(2,i),'r+','markersize',10.0,'linewidth',1.5);
        text(a(1,i)+delta,a(2,i)+delta,num2str(i),'fontsize',16,'fontweight','bold');
    elseif Ndim==3,
        plot3(a(1,i),a(2,i),a(3,i),'r+','markersize',10.0,'linewidth',1.5);
        text(a(1,i)+delta,a(2,i)+delta,a(3,i)+delta,num2str(i),'fontsize',16,'fontweight','bold');
    end;
end;

if Ndim==2,
    plot(a(1,:),a(2,:),'b-','linewidth',3.0); 
    axis([min_a(1)-0.5, max_a(1)+0.5, min_a(2)-0.5, max_a(2)+0.5]);
elseif Ndim==3,
    plot3(a(1,:),a(2,:),a(3,:),'b-','linewidth',3.0); 
    axis([min_a(1)-0.5, max_a(1)+0.5, min_a(2)-0.5, max_a(2)+0.5,min_a(3)-0.5, max_a(3)+0.5]);
end;

set(gcf,'color','w');
set(gcf, 'unit', 'normalized', 'position', [0.1,0.1,0.8,0.8]);  % 设置窗口全屏
set(gca,'position',[0 0 1 1]); % 设置绘图区域（坐标系）在窗口中的位置, position = [left, bottom, width, height]
axis off equal;
