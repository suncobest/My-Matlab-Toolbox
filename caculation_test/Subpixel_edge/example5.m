%% SUBPIXEL EDGES - EXAMPLE 5 -----------
% SUBPIXEL EDGE DETECTION IN A HIGH-NOISE SYNTHETIC IMAGE

%% syntethic ring
addpath('Synthetic');
imageSize = 35;
xCenter = imageSize/2;
yCenter = imageSize/2;
innerRadius = 8.0;
outerRadius = 10.0;
innerIntensity = 100;
outerIntensity = 200;
gridResolution = 100;
image = ring(imageSize, imageSize, xCenter, yCenter, ...
    innerRadius, outerRadius, innerIntensity, outerIntensity, ...
    gridResolution);

%% add noise
noisePercent = 5;
image = noise(image, noisePercent); 

%% subpixel detection
threshold = 20;
iter = 10;
I = image;

% regular form
% [edges, I] = subpixelEdges(I, threshold, 'SmoothingIter', iter);

%% show image
showRestoredImage = true;
if showRestoredImage
    imshow(I/255,'InitialMagnification', 'fit');
else
    imshow(image/255,'InitialMagnification', 'fit');
end

% step by step form
for n=1:iter
    [edges, I] = subpixelEdges(I, threshold);
    if showRestoredImage,
        imshow(I/255,'InitialMagnification', 'fit');
    else
        imshow(image/255,'InitialMagnification', 'fit');
    end;
    pause(0.5);
end

%% show edges
visEdges(edges);
