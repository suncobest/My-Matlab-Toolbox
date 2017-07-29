% show wings' kinematics displayed in one period
nx = 1280;
ny = 800;
height = 4;

Nb = length(ind_down)-1;
dfpc = round(ndiv/4);
XX = NaN( 6, ndiv+1, Nb);
id_down = NaN(2,Nb+1);
nm = (nn-1)/(N_active-1);
t0 = (0:ndiv)/ndiv;
id = round((ind_down(1)-ind_active(1))*nm)+1;
jj = max(1, id-dfpc);
ind = jj : id+dfpc;
[~,id] = min(strphi(:,ind),[],2);
id_down(:,1) = id+jj-1;
for ii=1:Nb,
    id = round((ind_down(ii+1)-ind_active(1))*nm)+1;
    jj = id-dfpc;
    ind = jj : min(id+dfpc,nn);
    [~,id] = min(strphi(:,ind),[],2);
    id_down(:,ii+1) = id+jj-1;
    ind = id_down(1,ii) : id_down(1,ii+1);
    id = (ind-ind(1))/(ind(end)-ind(1));
    XX([1,3,5], :, ii) = spline_interp3(id, [strphi(1,ind); strtheta(1,ind); strpsi(1,ind)], t0);
    ind = id_down(2,ii) : id_down(2,ii+1);
    id = (ind-ind(1))/(ind(end)-ind(1));
    XX([2,4,6], :, ii) = spline_interp3(id, [strphi(2,ind); strtheta(2,ind); strpsi(2,ind)], t0);
end;
a = mean(XX,3);
b = sqrt(sum((XX-a(:,:,ones(1,Nb))).^2, 3)/(Nb-1));
avphi = a(1:2,:);
avtheta = a(3:4,:);
avpsi = a(5:6,:);
stdphi = b(1:2,:);
stdtheta = b(3:4,:);
stdpsi = b(5:6,:);

[a, id] = max(XX(1:2,:,:),[],2);
maxphi = mean(a,3);
Tmaxphi = (mean(id,3)-1)/(ndiv+1);
minphi = mean(XX(1:2,1,:),3);

a = [0,1];
m = 3;
n = 2;
d = 0.1;    % margin of figure
w = (1-(n+0.5)*d)/n;    % width of figure
e = a(1)+[-0.4, 0.1]*d*(a(2)-a(1))/w;   % x postion of ylabel and text
h = (1-(m+2)*d/2)/m;    % height of figure
c = (0:n-1)*(w+d)+d;
d = (0:m-1)*(h+d/2)+d;
w1 = 0.5;
palette = lines(7);

t = [t0,t0(ndiv+1:-1:1)];
w2 = 0.2;       % alpha of fill function
m = 14;      % font size
figure;
% draw average wing kinematics
XX = [avpsi; avpsi+stdpsi; avpsi-stdpsi];
b = [min(XX(:)), max(XX(:))];
f = (b(2)-b(1))/10;
b = round(b+[-1,1]*f);
axes('position',[c(1) d(1) w h]);
plot(Tmaxphi(2)*[1,1], b, 'k-');
hold on;
temp = [XX(3:4,:), XX(5:6,ndiv+1:-1:1)];
plot(t0', XX(2,:)', '-', 'color', palette(3,:), 'linewidth',w1);
fill(t,temp(2,:),palette(3,:),'EdgeColor','none','FaceAlpha',w2);
xlabel('$\hat t$','Interpreter','latex','fontname','Arial','fontsize',m);
ylabel('$\overline{\psi_l}\ \ (^\circ)$','Interpreter','latex','fontname','Arial','fontsize',m,'position',[e(1), sum(b)/2]);
text(e(2), b(2)-f/2, ['({\it\fontname{Arial}\fontsize{' num2str(m) '}c})']);
set(gca,'YMinorTick','on','box','off','xlim',a,'ylim',b);

axes('position',[c(2) d(1) w h]);
plot(Tmaxphi(1)*[1,1], b, 'k-');
hold on;
plot(t0', XX(1,:)','-', 'color', palette(2,:), 'linewidth', w1);
fill(t,temp(1,:),palette(2,:),'EdgeColor','none','FaceAlpha',w2);
xlabel('$\hat t$','Interpreter','latex','fontname','Arial','fontsize',m);
ylabel('$\overline{\psi_r}\ \ (^\circ)$','Interpreter','latex','fontname','Arial','fontsize',m,'position',[e(1), sum(b)/2]);
text(e(2), b(2)-f/2, ['({\it\fontname{Arial}\fontsize{' num2str(m) '}f})']);
set(gca,'YMinorTick','on','box','off','xlim',a,'ylim',b);


XX = [avtheta; avtheta+stdtheta; avtheta-stdtheta];
b = [min(XX(:)), max(XX(:))];
f = (b(2)-b(1))/10;
b = round(b+[-1,1]*f);
axes('position',[c(1) d(2) w h]);
plot(Tmaxphi(2)*[1,1], b, 'k-');
hold on;
plot(t0', XX(2,:)', '-', 'color', palette(3,:), 'linewidth',w1);
temp = [XX(3:4,:), XX(5:6,ndiv+1:-1:1)];
fill(t,temp(2,:),palette(3,:),'EdgeColor','none','FaceAlpha',w2);
ylabel('$\overline{\theta_l}\ \ (^\circ)$','Interpreter','latex','fontname','Arial','fontsize',m,'position',[e(1), sum(b)/2]);
text(e(2), b(2)-f/2, ['({\it\fontname{Arial}\fontsize{' num2str(m) '}b})']);
set(gca,'XTickLabel',[],'xcolor','k','YMinorTick','on','box','off','xlim',a,'ylim',b);


axes('position',[c(2) d(2) w h]);
plot(Tmaxphi(1)*[1,1], b, 'k-');
hold on;
plot(t0', XX(1,:)','-', 'color', palette(2,:), 'linewidth', w1);
fill(t,temp(1,:),palette(2,:),'EdgeColor','none','FaceAlpha',w2);
ylabel('$\overline{\theta_r}\ \ (^\circ)$','Interpreter','latex','fontname','Arial','fontsize',m,'position',[e(1), sum(b)/2]);
text(e(2), b(2)-f/2, ['({\it\fontname{Arial}\fontsize{' num2str(m) '}e})']);
set(gca,'XTickLabel',[],'xcolor','k','YMinorTick','on','box','off','xlim',a,'ylim',b);


XX = [avphi; avphi+stdphi; avphi-stdphi];
b = [min(XX(:)), max(XX(:))];
f = (b(2)-b(1))/10;
b = round(b+[-1,1]*f);
axes('position',[c(1) d(3) w h]);
plot(Tmaxphi(2)*[1,1], b, 'k-');
hold on;
plot(t0', XX(2,:)', '-', 'color', palette(3,:), 'linewidth',w1);
temp = [XX(3:4,:), XX(5:6,ndiv+1:-1:1)];
fill(t,temp(2,:),palette(3,:),'EdgeColor','none','FaceAlpha',w2);
ylabel('$\overline{\phi_l}\ \ (^\circ)$','Interpreter','latex','fontname','Arial','fontsize',m,'position',[e(1), sum(b)/2]);
text(e(2), b(2)-f/2, ['({\it\fontname{Arial}\fontsize{' num2str(m) '}a})']);
set(gca,'XTickLabel',[],'xcolor','k','YMinorTick','on','box','off','xlim',a,'ylim',b);


axes('position',[c(2) d(3) w h]);
plot(Tmaxphi(1)*[1,1], b, 'k-');
hold on;
plot(t0', XX(1,:)','-', 'color', palette(2,:), 'linewidth', w1);
fill(t,temp(1,:),palette(2,:),'EdgeColor','none','FaceAlpha',w2);
ylabel('$\overline{\phi_r}\ \ (^\circ)$','Interpreter','latex','fontname','Arial','fontsize',m,'position',[e(1), sum(b)/2]);
text(e(2), b(2)-f/2, ['({\it\fontname{Arial}\fontsize{' num2str(m) '}d})']);
set(gca,'XTickLabel',[],'xcolor','k','YMinorTick','on','box','off','xlim',a,'ylim',b);

resolution = round(ny/height);   % dpi
set(gcf,'color','w','PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]*2/resolution);
drawnow;
save_name = [imgdir '/Av_eulerAngle'];
print(gcf,'-djpeg',['-r' num2str(resolution*2)],save_name);
saveas(gcf, save_name,'fig');
