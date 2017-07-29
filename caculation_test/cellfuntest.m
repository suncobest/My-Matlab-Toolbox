function result = cellfuntest
a={{1:5;3:10},{4:9;2:8}};
result = cellfun(@maxtimes,a);
end

function mt = maxtimes(input)
mt = max(input{1})*max(input{2});
end
