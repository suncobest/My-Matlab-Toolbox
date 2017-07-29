%% gaussian pyramid

img1=imread('baboon200.jpg');
[m,n]=size(img1);
w=fspecial('gaussian',[3 3]);
img2=imresize(imfilter(img1,w),[m/2 n/2]);
img3=imresize(imfilter(img2,w),[m/4 n/4]);
img4=imresize(imfilter(img3,w),[m/8 n/8]);
img5=imresize(imfilter(img4,w),[m/16 n/16]);
imshow(img1);
figure,imshow(img2);
figure,imshow(img3);
figure,imshow(img4);
figure,imshow(img5);