%% 函数句柄/function_handle(@)
%   句柄是一种间接调用函数的方式。
% 
% 语法
%    handle=@functionname
% 
%    handle=@(arglist)anonymous_function

min_and_max = @(x) cellfun(@(f) f(x), {@min, @max});
% extremum = min_and_max (randi(10,8,1))
[extremum,indices] = min_and_max (randi(10,8,1))

%% 映射
% 在上面，我们创建了一个关于x的映射，也就是一个函数。下面，我们可以创建一个映射函数
% 来对多个输入分别进行多种函数操作，听起来挺复杂，看懂了也就很直观了。我们定义val为
% 一个cell数组，里面包含了我们的输入数据，fcns是一系列函数句柄。

map = @(val, fcns) cellfun(@(f) f(val{:}), fcns); 

[[extremum, indices] = map({x}, {@min, @max});