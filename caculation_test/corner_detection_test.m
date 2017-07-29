imName = 'corner_gray_small.png';
root = 'D:\My Documents\zpf\文章\特征提取';
I = imread(fullfile(root,imName));
if ndims(I)==3
    I = 0.299 * I(:,:,1) + 0.587 * I(:,:,2) + 0.114 * I(:,:,3);  % 对于RGB图片只保留亮度信息（YUV的Y分量），对于灰度图像，则省去该步
end
%% cornerfinder 测试（基于Harris算法）
imshow(I),hold on
xt = ginput(6)';
[xc,good,bad,type] = cornerfinder(xt,I,7,7)
plot(xt(1,:),xt(2,:),'r+')
plot(xc(1,:),xc(2,:),'g.')

%% Harris 角点检测函数
% MinQuality表示角点的最低质量，取值[0,1]，越大角点质量越好，可以通过调大阈值删除错误角点；FilterSize表示高斯滤波器的窗口大小，领域必须为奇数，取值[3,inf)，
% 高斯滤波器的标准差为FilterSize/3；ROI表示兴趣区域[x,y,w,h]，前两位数表示兴趣区域的左上角点[x,y]，即行列标[j,i]；后两位表示矩形宽和高。 
% points = detectHarrisFeatures(I,'MinQuality',0.2,'FilterSize',5,'ROI', [1,1,size(I,2),size(I,1)]);
points = detectHarrisFeatures(I,'MinQuality',0.02,'FilterSize',5,'ROI', [50,150,100,100]);
strongest = points.selectStrongest(50);
imshow(I); hold on;
plot(points);
X = strongest.Location
points.Location
points.Metric
strongest.Metric
plot(X(:,1), X(:,2), 'r+');

%% 按照函数detectHarrisFeatures定义的一维高斯掩模
vecNorm = @(x) x/sum(x(:));
GaussFilter = @(n,sigma) vecNorm( exp(-(-(n-1)/2:(n-1)/2).^2/(2*sigma.^2)) );    % @(X,sigma) vecNorm( exp(-X.^2/(2*sigma.^2)) );

gaussint = @(t) integral(@(x) exp(-x.^2/2)/sqrt(2*pi),-t,t);


%% cornerPoints 创建对象,此对象拥有3个参数（Location，Count，Metric）和4个方法（isempty，plot，selectStrongest，size）
% Location表示角点坐标，Count表示个数，Metric表示角点强度的度量

I = checkerboard(50,2,2);
location = [51    51    51   100   100   100   151   151   151;...
            50   100   150    50   101   150    50   100   150]';
points = cornerPoints(location);
imshow(I); hold on;
points.plot  % plot(points);
points.Location
points.Count
points.Metric


%% Harris 角点检测函数
I=checkerboard(100,2,3);
points = detectHarrisFeatures(I);
strongest = points.selectStrongest(20);
imshow(I); hold on;
plot(points);

%% FAST 角点检测函数（角点与周围亮度有明显差别，且领域的点能够连起来超过3/4圆弧）
% MinContrast表示角点与其周围像素的最小亮度差，取值（0,1），增加此数值可减少探测的角点数
% points = detectFASTFeatures(I,'MinQuality',0.2,'MinContrast',0.3,'ROI', [30,60,500,300]);
points = detectFASTFeatures(I);
strongest = points.selectStrongest(54);
imshow(I); hold on;
plot(strongest);

%% detectMinEigenFeatures 最小特征值角点检测函数（Shi and Tomasi ）与Harris方法类似
% points = detectMinEigenFeatures(I);
points = detectMinEigenFeatures(I,'MinQuality',0.02,'FilterSize',5,'ROI', [50,150,100,100]);
strongest = points.selectStrongest(54);
imshow(I); hold on;
plot(strongest);

%% detectSURFFeatures函数
I = imread('cameraman.tif');
points = detectSURFFeatures(I);
imshow(I); hold on;
plot(points.selectStrongest(10));

%% detectMSERFeatures
I = imread('cameraman.tif');
regions = detectMSERFeatures(I);
imshow(I); hold on;
plot(regions);
