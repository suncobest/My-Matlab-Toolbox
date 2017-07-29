% 证明png无损

for kk = 1:10
    I = rand(720,1280,3);
    I = im2uint8(I);
    eval(['I' num2str(kk) ' = I;']);
    imwrite(I, ['I' num2str(kk) '.png'],'png')    
end

% 读取图片并与原变量比较
for kk=1:10
    I = imread(['I' num2str(kk) '.png']);
    eq = (I == eval(['I' num2str(kk)]));
    eq = find(eq==0);  % 找出不相等的元素
    eval(['eq' num2str(kk) ' = eq;']);
end