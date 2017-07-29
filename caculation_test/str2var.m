% 字符串转换成变量名

% eval：把字符串当命令来执行
%  例子：
clear,clc,echo off
NameSource=['a','b','c'];
for i=1:3
    Name(i)={['VarName',num2str(i)]}; % 创建cell
    eval([Name{i},'=NameSource(i)'])
end;

%% 变量名转换成字符串

% who返回变量名，例：
a1=123;a2=1234;a3=444;
 b=who('a*')    %返回的b是cell
