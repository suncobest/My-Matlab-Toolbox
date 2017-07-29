clear;
figure;
for i=1:25
[f,c]=square_wave(1,1,2*i);
subplot(5,5,i),plot(c,f,'b-'),title(['Order of ',num2str(2*i-1)]);
end
i=getframe(gcf);
colormap(i.colormap);
imwrite(i.cdata,'Square_wave.png')