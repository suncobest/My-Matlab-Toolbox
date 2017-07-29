imgdir='E:\InsectFlying\1506_5cam\selected\2_15061117\old';
imgbase = '3D_projection_';
imfile = dir([imgdir '/' imgbase '*.fig']);
assert(~isempty(imfile),'No fig file found with the specific base name in the directory!');
openfig([imgdir '/' imfile(1).name]);
flag = 1;
while flag,
    azel = input('Please set the aziumth and elevation of the 3D view: ([az, el]) ');
    if length(azel)~=2,
        warning('Unexpected input!');
        continue;
    end;
    figure(gcf);
    view(azel); 
    flag = input('Reset the orientation or not? ([]=yes, other=no) ','s');
    flag = isempty(flag);
end;
close(gcf);
for i=1:length(imfile),
    save_name = [imgdir '/' imfile(i).name(1:end-4)];
    open(save_name);
    view(azel);
    print(gcf,'-djpeg','-r300',save_name);
    saveas(gcf, save_name,'fig');
    close(gcf);
end;