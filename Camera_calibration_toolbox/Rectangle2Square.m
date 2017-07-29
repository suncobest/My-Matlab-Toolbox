function Square = Rectangle2Square(rectangle,L,W)

% Generate the square from a rectangle of known segment lengths
% from pt1 to pt2 : L
% from pt2 to pt3 : W

[u_hori,u_vert] = UnWarpPlane(rectangle);   % 由角点的像生成三维网格矢量，相当于将实际网格进行缩放，直到第一个角点缩放到归一化平面

% n_u_hori and n_u_vert are unit vectors,scale is a zooming factor
% pt1 = pt1;
% pt2 = pt1+u_hori = pt1+scale*L*n_u_hori;
% pt3 = pt1+u_hori+u_vert = pt1+scale*L*n_u_hori+scale*W*n_u_vert;
% pt4 = pt1+v_ert = pt1+scale*W*n_u_vert;


coeff_x = sqrt(W/L);  % 不改变面积，由长L宽W的长方形变成边长sqrt(L*W)的正方形，第一个角点pt1保持不动
coeff_y = 1/coeff_x;

x_coord = [ 0 coeff_x  coeff_x 0];      % [0  sqrt(W/L)  sqrt(W/L)  0]
y_coord = [ 0 0 coeff_y coeff_y];       % [0  0  sqrt(L/W)  sqrt(L/W)]


X = rectangle(:,1) * ones(1,4) + u_hori*x_coord + u_vert*y_coord;

% X1 = pt1;
% X2 = pt1+scale*sqrt(L*W)*n_u_hori;
% X3 = pt1+scale*sqrt(L*W)*(n_u_hori+n_u_vert);
% X4 = pt1+scale*sqrt(L*W)*n_u_vert;


Square = X ./ (ones(3,1)*X(3,:));
% 将生成的正方形网格投影到归一化平面上,得到正方形网格的像

return;





%% draw
% keyboard;

delta = max(abs(X(:)))/50;
figure(1);
plot3(X(1,:),X(3,:),-X(2,:),'r.',X(1,[1:4 1]),X(3,[1:4 1]),-X(2,[1:4 1]),'c-');
for i =1:4
text(X(1,i)+delta,(X(3,i)+delta),-X(2,i)+delta,num2str(i),'color','c');
end;
plot3(Square(1,:),Square(3,:),-Square(2,:),'r.',Square(1,[1:4 1]),Square(3,[1:4 1]),-Square(2,[1:4 1]),'r-');
for i =1:4
text(Square(1,i)+delta,(Square(3,i)+delta),-Square(2,i)+delta,num2str(i),'color','r');
end;
axis tight

xos=reshape([X;zeros(3,4);Square],3,[]);
figure(1),plot3(xos(1,:),xos(3,:),-xos(2,:),'color',0.8*[1 1 1]);

return;




%% test: uncomment the 1st key word 'return' (also in 'UnWarpPlane.m') and run the following code in the workspace
X0 =[
     0    10    10     0
     0     0    20    20
     0     0     0     0];
R=rodrigues(rand(3,1));
Y=R*X0+repmat([-5;-5;10],1,4);
y=Y./(ones(3,1)*Y(3,:));
square = Rectangle2Square(y,10,20);
