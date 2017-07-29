[rt,f,g]=twosin_comp(1,1/(4*pi),3/(4*pi),320,240);
f=mat2gray(f);
g=mat2gray(g);
figure,imshow(f);axis equal;
figure,imshow(g);axis equal;