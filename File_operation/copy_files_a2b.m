source = 'G:\SetupFiles\Fonts\AdobeFonts';
destination = 'G:\SetupFiles\Fonts\Adobe';
docs=dir(source);
count = length(docs);
tmp = double('K');
for i = 1:count,    
    files = docs(i).name;
    if double(files(1)) < tmp,
        continue;
    end;        
    files = [source, '\', files];
    copyfile(files,destination);
end;