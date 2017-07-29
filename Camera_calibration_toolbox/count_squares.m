function ns = count_squares(I,x1,y1,x2,y2,win)

[ny,nx] = size(I);

if ((x1-win <= 0) || (x1+win >= nx) || (y1-win <= 0) || (y1+win >= ny) || ...
        (x2-win <= 0) || (x2+win >= nx) || (y2-win <= 0) || (y2+win >= ny))
    ns = -1;
    return;
end;

if ((x1 - x2)^2+(y1-y2)^2) <  win,
    ns = -1;
    return;
end;

% 此函数需要在lambda两侧不对称,且使得out_ff正负交替（很重要），才能数格子
% lambda is the line pass through homogeneous points [x1;y1;1] and [x2;y2;1]
lambda = [y1 - y2;x2 - x1;x1*y2 - x2*y1];                 % lambda = cross([x1;y1;1],[x2;y2;1])
lambda = 1/sqrt(lambda(1)^2 + lambda(2)^2) * lambda;      % 法式直线方程[a b c]*x/sqrt(a^2+b^2)

% l1 = lambda + [0;0;win];
% l2 = lambda - [0;0;win]; 
% 直线l1和l2平行于直线lambda，分别在其两侧，距离都是win。兴趣区域为l1和l2之间的带状区域

dx = x2 - x1;
dy = y2 - y1;

if abs(dx) > abs(dy),                                     % 直线lambda与x轴夹角小于45度，取y=f(x)；否则取x=f(y)
   if x2 > x1,
      xs = x1:x2;
   else
      xs = x1:-1:x2;
   end;
   ys = -(lambda(3) + lambda(1)*xs)/lambda(2);             % 点列[xs;ys;1]在直线lambda上
else
   if y2 > y1,
       ys = y1:y2;
   else
       ys = y1:-1:y2;
   end;
   xs = -(lambda(3) + lambda(2)*ys)/lambda(1);
end;

Np = length(xs);
xs_mat = ones(2*win + 1,1)*xs;
ys_mat = ones(2*win + 1,1)*ys;
win_mat = (-win:win)'*ones(1,Np);
xs_mat2 = round(xs_mat - win_mat * lambda(1));              % 构造点阵[xs,ys]-n*[a,b]，n=-win:win，其中norm([a,b])=1
ys_mat2 = round(ys_mat - win_mat * lambda(2));              % 得到与直线lambda距离-win到win的像素点坐标
ind_mat = (xs_mat2 - 1) * ny + ys_mat2;                     % 从像素坐标[x,y]转化为单个下标(x-1)*ny+y
ima_patch = zeros(2*win + 1,Np);
ima_patch(:) = I(ind_mat(:));                               % 带状区域的像素灰度

%ima2 = ima_patch(:,win+1:end-win);

filtk = [ones(win,Np);zeros(1,Np);-ones(win,Np)];           % 生成带状系数，上半部分（高win）为1，中间一行为0，下半部分为-1
out_f = sum(filtk.*ima_patch);                              % 将带状区域的像素灰度上半部分减去下半部分，再求和（沿垂直于lambda方向）
out_ff = conv2(out_f,[1/4 1/2 1/4],'same');                % out_f为向量，对一维向量卷积使用conv和conv2没有区别，线性帐篷tent模糊
out_ff = out_ff(win+1:end-win);                           % 去掉两个端点：在两头分别截去半径win的长度
ns = length(find(((out_ff(2:end)>=0)&(out_ff(1:end-1)<0)) | ((out_ff(2:end)<=0)&(out_ff(1:end-1)>0))))+1; 

% 寻找out_ff改变符号的点的个数（中间的角点数），然后+1得到格子数。
% (out_ff(2:end)>=0)&(out_ff(1:end-1)<0)找出n-1为负，n为非负的转折点（包括从-到+）
% (out_ff(2:end)<=0)&(out_ff(1:end-1)>0)找出n-1为正，n为非正的转折点（包括从+到-）
% 此函数不能找出从0到-的转折点，也不能找出从0到+的转折点

return;
