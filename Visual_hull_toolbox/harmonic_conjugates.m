%%  cross ratio (ABCD)=-1
normalize = @(x) x./(ones(size(x,1),1)*x(end,:)); 
combpoints = @(P1,P2,lamda1,lamda2) lamda1*P1+lamda2*P2;

% 调和共轭点对
ncp = @(P1,P2,lamda) normalize([combpoints(P1,P2,1,lamda),combpoints(P1,P2,1,-lamda)]);

% line是直线的齐次坐标，x为横坐标，y为纵坐标，nyfx由x输出y，nxfy由y输出x
% 建议当直线斜率的绝对值<1时，使用nxyFx；否则使用nxyFy
nxyFx = @(line,x) [x; -(x*line(1)+line(3))/line(2)];  
nxyFy = @(line,y) [-(y*line(2)+line(3))/line(1); y];  

A=[0;0;1];
B=[2;0;1];
E=[1.5;1.5;1];

AB=[0;1;0];
AE=cross(A,E);   % AE is the line passing the two points A and E. 
BE=cross(B,E);
F=[nxyFx(BE,1.75);1];  % F is a fixed point on BE
AF=cross(A,F);

tri=[A,B,E,A,F];

% extending lines
limit=20;
PAB = nxyFx(AB,limit*[-1,1]);   % 左右两个极限点的横坐标为limit*[-1,1]
PAE = nxyFx(AE,limit*[-1,1]);
PBE = nxyFx(BE,limit*[-1,1]);
PAF = nxyFx(AF,limit*[-1,1]);

%%
figure(1);
set(gcf,'Color','w');   % set(gcf,'Color',[1 1 1]);

xylimit=[-3,5,-1,3];

k1=E(2)/(E(1)-xylimit(1));
k2=E(2)/(E(1)-xylimit(2));

dt=0.1;
filename='Harmonic conjugates.gif';

for theta=0:2:179
    if theta == 90;
        EC=[1;0;-E(1)];       
    else
        EC=[-tand(theta); 1; E(1)*tand(theta)-E(2)];  % rad: tan; degree: tand
    end
    C=cross(EC,AB);
    G=cross(EC,AF);
    BG=cross(B,G);
    H=cross(BG,AE);
    FH=cross(F,H);
    D=cross(FH,AB);
    
    % extending lines (normalized points)
    if abs(EC(1)) > abs(EC(2))
        PEC = nxyFy(EC,limit*[-1,1]);
    else
        PEC = nxyFx(EC,limit*[-1,1]);
    end
    
    if abs(BG(1)) > abs(BG(2))
        PBG = nxyFy(BG,limit*[-1,1]);
    else
        PBG = nxyFx(BG,limit*[-1,1]);
    end
    
    if abs(FH(1)) > abs(FH(2))
        PFH = nxyFy(FH,limit*[-1,1]);
    else
        PFH = nxyFx(FH,limit*[-1,1]);
    end
    
    C=C/C(3);
    D=D/D(3);
    G=G/G(3);
    H=H/H(3);

%     C(1)
%     D(1)
%     Crossratio=C(1)*(D(1)-B(1))/((C(1)-B(1))*D(1))
    
    plot(PAB(1,:),PAB(2,:),PAE(1,:),PAE(2,:),PBE(1,:),PBE(2,:),PAF(1,:),PAF(2,:),'k-','LineWidth',1)
    hold on
    plot(tri(1,:),tri(2,:),'k-','LineWidth',2)
    plot(PEC(1,:),PEC(2,:),'b-',PBG(1,:),PBG(2,:),'r-',PFH(1,:),PFH(2,:),'g-','LineWidth',1.5)
    plot(C(1),C(2),'bo',D(1),D(2),'go','Markersize',5,'linewidth',2)
    plot(G(1),G(2),'ko',H(1),H(2),'ko','Markersize',4,'linewidth',1)
    
    axis equal off
    axis(xylimit);
    set(gca,'position',[0.1,0.1,0.8,0.8])
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [ind,map] = rgb2ind(im,256);
    if  theta==0
        imwrite(ind,map,filename,'gif','LoopCount',Inf,'DelayTime',dt);
    else
        imwrite(ind,map,filename,'gif','WriteMode','append','DelayTime',dt);
    end
    
%     pause
    %    print(1,'-dpng',['unit circle ' sprintf('%03d',i) '.png']);
    hold off
    
end

%%
theta = 60;
EC=[-tand(theta); 1; E(1)*tand(theta)-E(2)];  % rad: tan; degree: tand
C=cross(EC,AB);
G=cross(EC,AF);
BG=cross(B,G);
H=cross(BG,AE);
FH=cross(F,H);
I=cross(FH,EC);
D=cross(FH,AB);

% normalization
C=C/C(3);
G=G/G(3);
H=H/H(3);
I=I/I(3);
D=D/D(3);

figure(1);
plot(tri(1,:),tri(2,:),'k-','LineWidth',2);
hold on;
plot([xylimit(1),xylimit(2)],[0,0],'k-','LineWidth',1);

plot([E(1),C(1),G(1),I(1)],[E(2),C(2),G(2),I(2)],'b-','LineWidth',2);
plot([B(1),H(1),G(1)],[B(2),H(2),G(2)],'r-','LineWidth',2);
plot([F(1),G(1)],[F(2),G(2)],'k-','LineWidth',1);
plot([A(1),H(1)],[A(2),H(2)],'k-','LineWidth',1);
plot([F(1),D(1),H(1),I(1)],[F(2),D(2),H(2),I(2)],'g-','LineWidth',2);

plot(C(1),C(2),'ro','Markersize',6,'linewidth',1.5);
plot(D(1),D(2),'ro','Markersize',6,'linewidth',1.5);

axis equal off;
axis(xylimit);
set(gca,'position',[0.1,0.1,0.8,0.8]);
set(gcf,'Color','w'); 

drawnow
frame = getframe(1);
im = frame2im(frame);
[ind,map] = rgb2ind(im,256);
imwrite(ind,map,'harmnonic conjugates.png','png')
 
% imwrite(background,[pathname 'background.png'],'png','alpha',alpha_trans)
% imwrite(repmat(background,1,1,3),gray(256),[pathname 'background.png'],'alpha',alpha_trans)
% print(1,'-dpng',['harmnonic conjugates.png'])

hold off

    

%%

% figure(1);
% set(gcf,'Color','w');   % set(gcf,'Color',[1 1 1]);
% set(gca,'position',[0.1,0.1,0.8,0.8]);
% 
% xylimit=[-2,4,-0.5,2];
% 
% k1=E(2)/(E(1)-xylimit(1));
% k2=E(2)/(E(1)-xylimit(2));
% 
% dt=0.5;
% filename='Harmonic conjugates.gif';

% for i=0:20
%     C=[i-10;0;1];
%     EC=cross(E,C);
%     G=cross(EC,AF);
%     BG=cross(B,G);
%     H=cross(BG,AE);
%     FH=cross(F,H);
%     D=cross(FH,AB);
%     I=cross(FH,EC);
%     G=G/G(3);
%     H=H/H(3);
%     D=D/D(3);
%     I=I/I(3);
%     
%     C(1)
%     D(1)
%     Crossratio=C(1)*(D(1)-B(1))/((C(1)-B(1))*D(1))
%     
%     plot(tri(1,:),tri(2,:),'k-','LineWidth',2);
%     hold on;
%     
%     plot([A(1),C(1),D(1)],[A(2),C(2),D(2)],'k-','LineWidth',1);
%     plot([E(1),C(1),G(1),I(1)],[E(2),C(2),G(2),I(2)],'b-','LineWidth',2);
%     plot([B(1),H(1),G(1)],[B(2),H(2),G(2)],'r-','LineWidth',2);
%     plot([F(1),G(1)],[F(2),G(2)],'k-','LineWidth',1);
%     plot([A(1),H(1)],[A(2),H(2)],'k-','LineWidth',1);
%     plot([F(1),D(1),H(1),I(1)],[F(2),D(2),H(2),I(2)],'g-','LineWidth',2);
%     
%     plot(C(1),C(2),'ro','Markersize',6,'linewidth',1.5);
%     plot(D(1),D(2),'ro','Markersize',6,'linewidth',1.5);
%     
%     axis equal off;
%     drawnow
%     frame = getframe(1);
%     im = frame2im(frame);
%     [ind,map] = rgb2ind(im,256);
%     if  ~i
%         imwrite(ind,map,filename,'gif','LoopCount',Inf,'DelayTime',dt);
%     else
%         imwrite(ind,map,filename,'gif','WriteMode','append','DelayTime',dt);
%     end
%     
%     pause
%     hold off
%     
% end