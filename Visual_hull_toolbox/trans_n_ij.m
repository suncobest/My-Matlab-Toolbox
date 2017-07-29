function l=trans_n_ij(m)

% l=trans_n_ij(m) transform array (i,j) to the corresponding number n; or
% number n to the corresponding array (i,j). There is a one-to-one
% correspondence between array (i,j) and number n.

switch length(m)
    case 1
        l=n2ij(m);
    case 2
        l=ij2n(m);
    otherwise
        disp('Unexpected input! The input must either be a 2D array or a number!');
        return;
end


function l=n2ij(m)

r=round(sqrt(2*m));
if mod(r,2)
    j=m-r*(r-1)/2;
    i=r+1-j;
else
    i=m-r*(r-1)/2;
    j=r+1-i;
end

l=[i,j];
return


function l=ij2n(m)

i=m(1);
j=m(2);
r=i+j-1;
if mod(r,2)
    l=j+r*(r-1)/2;   
else
    l=i+r*(r-1)/2; 
end
return
