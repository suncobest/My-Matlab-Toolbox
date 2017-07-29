%% SUBPIXEL EDGES - EXAMPLE 6 -----------
% SUBPIXEL EDGE DETECTION IN A REAL ANGIOGRAPHY

%% load image
% url='http://serdis.dis.ulpgc.es/~atrujill/ngImgCrop-master/test/angio2.PNG';
url='angio2.PNG';
image = rgb2gray(imread(url));

%% subpixel detection
threshold = 4;
iter = 3;
[edges, RI] = subpixelEdges(image, threshold, 'SmoothingIter', iter); 

%% show image
showRestoredImage = false;
if showRestoredImage
    imshow(I/255,'InitialMagnification', 'fit');
else
    imshow(image,'InitialMagnification', 'fit');
end

%% show edges
visEdges(edges);

