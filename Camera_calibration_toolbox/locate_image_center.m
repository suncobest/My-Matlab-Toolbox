function cc_out = locate_image_center(cc_in, nxy1, nxy2 )
%locate_image_center relocate image center if image have different size.
assert(isvector(cc_in) && length(cc_in)==2, 'Unexpected input for the 1st argument!');
assert(isvector(nxy1) && length(nxy1)==2, 'Unexpected input for the 2nd argument!');
cc_out = cc_in;
if nargin<3,  % do not resize
    return;
else
    assert(isvector(nxy2) && length(nxy2)==2, 'Unexpected input for the 3rd argument!');
end;
nxy1 = nxy1(:);
nxy2 = nxy2(:);
if isequal(nxy1, nxy2),
    return;
end;
cc_out(:) = cc_out(:) + (nxy2-nxy1)/2;

return;

