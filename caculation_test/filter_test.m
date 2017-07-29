I = imread('cameraman.tif');
subplot(2,2,1); 
imshow(I); title('Original Image');

H = fspecial('motion',20,45);
MotionBlur = imfilter(I,H,'replicate');
subplot(2,2,2); 
imshow(MotionBlur);title('Motion Blurred Image');

H = fspecial('disk',10);
blurred = imfilter(I,H,'replicate');
subplot(2,2,3); 
imshow(blurred); title('Blurred Image');

H = [0 -1 0;-1 5 -1;0 -1 0];
blurred = imfilter(I,H,'replicate');
subplot(2,2,4);
imshow(blurred); title('Sharppened Image');

%%
clc
i = rand(5,5);
h = 1:5;
filter2(h,i,'full')
% conv(i,fliplr(h))
conv2(i,fliplr(h))
xcorr2(i,h)


