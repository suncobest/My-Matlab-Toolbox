% ֤��png����

for kk = 1:10
    I = rand(720,1280,3);
    I = im2uint8(I);
    eval(['I' num2str(kk) ' = I;']);
    imwrite(I, ['I' num2str(kk) '.png'],'png')    
end

% ��ȡͼƬ����ԭ�����Ƚ�
for kk=1:10
    I = imread(['I' num2str(kk) '.png']);
    eq = (I == eval(['I' num2str(kk)]));
    eq = find(eq==0);  % �ҳ�����ȵ�Ԫ��
    eval(['eq' num2str(kk) ' = eq;']);
end