% find the codedir of the toolbox
codeDir = mfilename('fullpath');
flag = find(codeDir=='\' | codeDir=='/',1,'last');
if ~isempty(flag),
    codeDir = codeDir(1 : flag(end)-1);
end;
if exist([codeDir '/load_imgdir.m'],'file'),
    load_imgdir;
    fprintf(1,'\nModify file ''load_imgdir.m'' in code directory if you want to redirect image path!\n');
else
    fprintf(1,'\nNo file ''load_imgdir.m'' found in the code directory!\nPlease compute visual hull first!\n');
    return;
end;

save_name = [imgdir '/kinematic_animation.mat'];
if exist(save_name,'file')==2,
    fprintf(1,'\nLoading kinematics parameters file ''kinematic_animation.mat'' ...\n');
    load(save_name);
else
    fprintf(1,'\nERROR: ''kinematic_animation.mat'' not found!\n');
    return;
end;

a = [ts(1),ts(end)]*1e3;
n = 6;
d = 0.1;    % margin of figure
w = 1-d*1.5;    % width and height of figure
h = (w-(n-1)*d/2)/n;
c = h+d/2;
c = (0:5)*c+d;
w1 = 1;
w2 = 0.5;
m = 9;      % font size
e = a(1)-(a(2)-a(1))*d/(2*w);   % x postion of ylabel

figure;
% draw wing orientation kinematics
XX = [strpsi; strpsi1];
b = [min(XX(:)), max(XX(:))];
b = round(b+[-1,1]*(b(2)-b(1))/10);
axes('position',[d c(1) w h]);
% plot(1e3*ts', XX(4,:)', '-', 'color', palette(7,:), 'linewidth',w2);
plot(1e3*ts', XX(4,:)', 'k-', 'linewidth',w2);
hold on;
plot(1e3*ts', XX(2,:)', '-', 'color', palette(3,:), 'linewidth', w1);
% plot(1e3*ts', XX(4,:)','k-',1e3*ts', XX(2,:)','r-');
xlabel('\fontsize{12}time (ms)');
% ylabel('\fontsize{12}{\it\psi_l} (^\circ)','position',[e, sum(b)/2]);
ylabel('$\psi_l\ \ (^\circ)$','Interpreter','latex','fontname','Arial','fontsize',m,'position',[e, sum(b)/2]);
set(gca,'YMinorTick','on','box','off','xlim',a,'ylim',b);
% legend('m','p');

axes('position',[d c(2) w h]);
% plot(1e3*ts', XX(3,:)','-', 'color', palette(7,:), 'linewidth', w2);
plot(1e3*ts', XX(3,:)','k-', 'linewidth',w2);
hold on;
plot(1e3*ts', XX(1,:)','-', 'color', palette(2,:), 'linewidth', w1);
% plot(1e3*ts', XX(3,:)','k-',1e3*ts', XX(1,:)','g-');
% ylabel('\fontsize{12}{\it\psi_r} (^\circ)','position',[e, sum(b)/2]);
ylabel('$\psi_r\ \ (^\circ)$','Interpreter','latex','fontname','Arial','fontsize',m,'position',[e, sum(b)/2]);
set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
% legend('m','p');

XX = [strtheta; strtheta1];
b = [min(XX(:)), max(XX(:))];
b = round(b+[-1,1]*(b(2)-b(1))/10);
axes('position',[d c(3) w h]);
% plot(1e3*ts', XX(4,:)', '-', 'color', palette(7,:), 'linewidth',w2);
plot(1e3*ts', XX(4,:)', 'k-', 'linewidth',w2);
hold on;
plot(1e3*ts', XX(2,:)', '-', 'color', palette(3,:), 'linewidth',w1); 
% plot(1e3*ts', XX(4,:)','k-',1e3*ts', XX(2,:)','r-');
% ylabel('\fontsize{12}{\it\theta_l} (^\circ)','position',[e, sum(b)/2]);
ylabel('$\theta_l\ \ (^\circ)$','Interpreter','latex','fontname','Arial','fontsize',m,'position',[e, sum(b)/2]);
set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
% legend('m','p');

axes('position',[d c(4) w h]);
% plot(1e3*ts', XX(3,:)', '-', 'color', palette(7,:), 'linewidth',w2);
plot(1e3*ts', XX(3,:)','k-', 'linewidth',w2);
hold on;
plot(1e3*ts', XX(1,:)','-', 'color', palette(2,:), 'linewidth',w1); 
% plot(1e3*ts', XX(3,:)','k-',1e3*ts', XX(1,:)','g-');
% ylabel('\fontsize{12}{\it\theta_r} (^\circ)','position',[e, sum(b)/2]);
ylabel('$\theta_r\ \ (^\circ)$','Interpreter','latex','fontname','Arial','fontsize',m,'position',[e, sum(b)/2]);
set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
% legend('m','p');

XX = [strphi; strphi1];
b = [min(XX(:)), max(XX(:))];
b = round(b+[-1,1]*(b(2)-b(1))/10);
axes('position',[d c(5) w h]);
% plot(1e3*ts', XX(4,:)','-', 'color', palette(7,:), 'linewidth',w2);
plot(1e3*ts', XX(4,:)', 'k-', 'linewidth',w2);
hold on;
plot(1e3*ts', XX(2,:)','-', 'color', palette(3,:), 'linewidth',w1);
% plot(1e3*ts', XX(4,:)','k-',1e3*ts', XX(2,:)','r-');
% ylabel('\fontsize{12}{\it\phi_l} (^\circ)','position',[e, sum(b)/2]);
ylabel('$\phi_l\ \ (^\circ)$','Interpreter','latex','fontname','Arial','fontsize',m,'position',[e, sum(b)/2]);
set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
% legend('m','p');

axes('position',[d c(6) w h]);
% plot(1e3*ts', XX(3,:)','-', 'color', palette(7,:), 'linewidth',w2);
plot(1e3*ts', XX(3,:)','k-', 'linewidth',w2);
hold on;
plot(1e3*ts', XX(1,:)','-', 'color', palette(2,:), 'linewidth',w1);
% plot(1e3*ts', XX(3,:)','k-',1e3*ts', XX(1,:)','g-');
% ylabel('\fontsize{12}{\it\phi_r} (^\circ)','position',[e, sum(b)/2]);
ylabel('$\phi_r\ \ (^\circ)$','Interpreter','latex','fontname','Arial','fontsize',m,'position',[e, sum(b)/2]);
set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
% legend('m','p');

% time of reversal
[~,idk] = min(strphi(1,1:ndiv));
Treversal = ts(idk) : 1/(2*cyclePS) : ts(end);
tupend = Treversal(3:2:end);    % end time of upstoke
Nb = length(tupend);
b = [tupend-1/(2*cyclePS); tupend]*1e3;	% period of upstoke
% plot transparent color bar on upstrokes
axes('position',[d, d, w, w]); hold on;
for ii=1:Nb,
    c = [b(1,ii), b(2,ii), b(2,ii), b(1,ii); 0, 0, 1, 1];
    fill(c(1,:),c(2,:),[1 1 1]*0.2,'EdgeColor','none','FaceAlpha',0.3);
end;
set(gca,'xlim',a); axis off;
resolution = round(ny/height);   % dpi
set(gcf,'color','w','PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]*2/resolution);
drawnow;
save_name = [imgdir '/wings_eulerAngle'];
print(gcf,'-djpeg',['-r' num2str(resolution*2)],save_name);
saveas(gcf, save_name,'fig');

fprintf(1,'done...\n');
