function  [u_hori,u_vert] = UnWarpPlane(x1,x2,x3,x4)

% Recovers the two 3D directions of the rectangular patch x1x2x3x4
% x1 is the origin point, ie any point of planar coordinate (x,y) on the
% rectangular patch will be projected on the image plane at:
% x1 + x * u_hori + y * u_vert
%
% Note: u_hori and u_vert are also the two vanishing points.

% x1,x2,x3,x4Ϊ���������ϵ�еĵ㣬�����һ��ƽ���ϵĵ�
if nargin < 4,
   
   x4 = x1(:,4);
   x3 = x1(:,3);
   x2 = x1(:,2);
   x1 = x1(:,1);
   
end;


% Image Projection:
L1 = cross(x1,x2);
L2 = cross(x4,x3);
L3 = cross(x2,x3);
L4 = cross(x1,x4);

% Vanishing point:
V1 = cross(L1,L2);
V2 = cross(L3,L4);

% Horizon line: ideal line, ������ƽ����3D��ʵ����ƽ���ƽ������ƽ��Ľ���
H = cross(V1,V2);

if H(3) < 0, H  = -H; end;


H = H / norm(H);   % ��ʧ�ߵĵ�λ���ꡣ������ĺ���ʧ�ߵ�ƽ��ΪPov����HΪPov�ĵ�λ��ʸ��


X1 = x1 / dot(H,x1);    % d1 = dot(H,x1)Ϊx1��ƽ��Pov�ľ��룬 X1 = x1/d1ʹ��X1��Pov�ľ���Ϊ1;
X2 = x2 / dot(H,x2);    % d2 = dot(H,x2)Ϊx2��ƽ��Pov�ľ��룬 X2 = x2/d2ʹ��X2��Pov�ľ���Ϊ1;
X3 = x3 / dot(H,x3);    % d3 = dot(H,x3)Ϊx3��ƽ��Pov�ľ��룬 X3 = x3/d3ʹ��X3��Pov�ľ���Ϊ1;
X4 = x4 / dot(H,x4);    % d4 = dot(H,x4)Ϊx4��ƽ��Pov�ľ��룬 X4 = x4/d4ʹ��X4��Pov�ľ���Ϊ1;

% ��Ȼx1,x2,x3,x4��X1,X2,X3,X4͸�Ӷ�Ӧ����X1,X2,X3,X4����ƽ��ƽ����Povƽ�棨ƽ��������ƽ�棩���Ҿ���Ϊ1��

scale = X1(3);     %   scale = x13/d1

X1 = X1/scale;      % �������ź�X1(3)=1,��X1�ڹ�һ��ƽ���ϣ���x1,x2,x3,x4���ڹ�һ��ƽ���ϣ���X1��x1�غϡ�
X2 = X2/scale;      % X1,X2,X3,X4����ƽ���X1��ƽ����Povƽ�棬���X1,X2,X3,X4��������������
X3 = X3/scale;      % һ���ı���x1,x2,x3,x4�ָ�Ϊƽ���ı���X1,X2,X3,X4
X4 = X4/scale;      

u_hori = X2 - X1;  % ��������Ϊ����ʸ��
u_vert = X4 - X1;

return;








 

%% draw 
% keyboard;
x=[x1 x2 x3 x4];
X=[X1 X2 X3 X4];
delta = max(abs(x(:)))/50;
V1=V1/V1(3);
V2=V2/V2(3);
lv1=[x1 x2 V1 x3 x4];
lv2 =[x2 x3 V2 x1 x4];

figure(1),plot3(lv1(1,:),lv1(3,:),-lv1(2,:),'r.',lv1(1,:),lv1(3,:),-lv1(2,:),'color',0.8*[1 1 1]);
axis equal, hold on
figure(1),plot3(lv2(1,:),lv2(3,:),-lv2(2,:),'r.',lv2(1,:),lv2(3,:),-lv2(2,:),'color',0.8*[1 1 1]);
figure(1),plot3([V1(1,:) V2(1,:)],[V1(3,:) V2(3,:)],-[V1(2,:) V2(2,:)],'color',0.8*[1 1 1]);
text(V1(1)+delta,V1(3)+delta,-(V1(2)+delta),'V_1','color','k');
text(V2(1)+delta,V2(3)+delta,-(V2(2)+delta),'V_2','color','k');

ox=reshape([x;zeros(3,4);x],3,[]);
figure(1),plot3(0,0,0,'r+',ox(1,:),ox(3,:),-ox(2,:),'color',0.8*[1 1 1]);
xX=reshape([x;X;x],3,[]);
figure(1),plot3(xX(1,:),xX(3,:),-xX(2,:),'color',0.8*[1 1 1]);
vo=[V1 zeros(3,1) V2];
figure(1),plot3(vo(1,:),vo(3,:),-vo(2,:),'color',0.8*[1 1 1]);
text(delta,delta,-delta,'O','color','r');

figure(1),plot3(x(1,:),x(3,:),-x(2,:),'r+',x(1,[1:4 1]),x(3,[1:4 1]),-x(2,[1:4 1]),'b-','linewidth',2.0);
for i =1:4
    text(x(1,i)+delta,(x(3,i)+delta),-x(2,i)+delta,['x_',num2str(i)],'color','b');
end
figure(1),plot3(X(1,:),X(3,:),-X(2,:),'r.',X(1,[1:4 1]),X(3,[1:4 1]),-X(2,[1:4 1]),'g-','linewidth',2.0);
for i =1:4
    text(X(1,i)+delta,(X(3,i)+delta),-X(2,i)+delta,['X_',num2str(i)],'color','g');
end

axis tight;
cross(H,cross(u_hori,u_vert))


return;



%% test : uncomment the 1st key word 'return' and run the following code in the workspace 

y =[
   -1.0511   -0.7144    0.2222    0.4772
    0.6964   -0.1180   -0.4444    0.6795
    1.0000    1.0000    1.0000    1.0000];
[du,dv]=UnWarpPlane(y);