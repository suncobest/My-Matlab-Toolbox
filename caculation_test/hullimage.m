% d=logical(randi(2,1,10)-1) 生成逻辑向量


a=randi(100,20,2);
b=convex_hull(a')';
hold off; plot(a(:,1),a(:,2),'k.');axis tight square;axis([0 100 0 100])
hold on; plot(b(:,1),b(:,2),'r+')
X=zeros(101,101);
[x,y]=ind2sub(size(X),1:numel(X));
x=x-1;y=y-1;
X=inpolygon(x',y',b(:,1),b(:,2));
figure,imshow(rot90(reshape(X,101,101)))  % 将矩阵逆时针旋转90°，从而将ij像素坐标系转换到常用xy坐标系
