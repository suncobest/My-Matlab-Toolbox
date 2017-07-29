%% SUBPIXEL EDGES - EXAMPLE 3 -----------
% SUBPIXEL EDGE DETECTION IN A BIG IMAGE OF A PRINTED TEXT CAPTURED BY A
% MOBILE PHONE CAMERA

%% load image
disp('Reading big image (8Mpixel)...');
% url='http://serdis.dis.ulpgc.es/~atrujill/subpixel/Portada%20paper.jpg';
url='Portada paper.jpg';
image = rgb2gray(imread(url));
imshow(image,'InitialMagnification', 'fit');

%% subpixel detection
disp('Computing subpixel edges...');
threshold = 10;
edges = subpixelEdges(image, threshold); 

%% show edges
visEdges(edges, 'showNormals', false);
disp('Done!');
