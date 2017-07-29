function u = demo_acwe(Img, iterNum)
%Effect: demonstrate the active contour without edge algorithm
%Inputs:
%Img: the input image
%iterNum: number of iteration
%Outputs:
%u: the final level-set function phi
%Author: Su dongcai at 2012/1/12 Email: suntree4152@gmail.com, qq:272973536
Img = double(Img(:, :, 1));
[ny,nx] = size(Img);
%apply median filter to denoise
Img = medfilt2(Img, [5, 5]);
%setting the initial level set function 'u':
center_len = [0.1,0.05];
center_len = round(center_len.*[ny,nx]);
lucorner = round([ny,nx].*[0.6,0.4]);
c0=2;
u = ones(ny, nx)*c0;
u(lucorner(1):lucorner(1)+center_len(1), lucorner(2):lucorner(2)+center_len(2))=-c0; 
%setting the parameters in ACWE algorithm:
mu=1;
lambda1=1; lambda2=1;
timestep = .1; v=1; epsilon=1;
%show the initial 0-level-set contour:
figure;imshow(Img, []);hold on;axis off,axis equal
title('Initial contour');
[c,h] = contour(u,[0 0],'r');
pause(0.1);
% start level set evolution
for n=1:iterNum
    u=acwe(u, Img,  timestep,...
             mu, v, lambda1, lambda2, 1, epsilon, 1);
    if mod(n,10)==0
        pause(0.1);
        imshow(Img, []);hold on;axis off,axis equal
        [c,h] = contour(u,[0 0],'r');
        iterNum=[num2str(n), ' iterations'];
        title(iterNum);
        hold off;
    end
end
imshow(Img, []);hold on;axis off,axis equal
[c,h] = contour(u,[0 0],'r');
totalIterNum=[num2str(n), ' iterations'];
title(['Final contour, ', totalIterNum]);

figure;
imagesc(u);axis off,axis equal;
title('Final level set function');