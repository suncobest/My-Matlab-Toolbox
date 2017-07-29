function y = sumn(n)
if n == 1
    y = 1 ;
else
    y = n + sumn(n-1) ;
end