function RenameSequence(old_base, file_ext, new_base, filepath, numfmt, keep_spacing)
% Rename file sequence: old_base to new_base, number to the format specified by numfmt.
% old_base: basename of file without number nor suffix.
% file_ext: suffix of files, containing all characters after numbers.
% new_base: new basename for files. (Default: keep the old basename).
% filepath: path of files to be renamed.
% numfmt: format of number: [start, ndigit] to set start and digit of number.
% keep the number spacing between file sequence or not£¬1 yes (default)£¬0 no (then spacing is 1).

%--------------------------------------------------------------------------
if nargin<6 || isempty(keep_spacing),
    keep_spacing = 1;
else
    keep_spacing = ~~keep_spacing ;
end;

if nargin<5 || isempty(numfmt),
    numfmt=[1, 4];      % [start, ndigit]
else
    numfmt = round(numfmt);
end;

if nargin<4 || isempty(filepath),
    filepath = pwd;
end;

if nargin<3 || isempty(new_base),
    new_base = old_base;
end;

if nargin<2,
    file_ext = [];      % most likely directory
end;

%filepath initialization   
if filepath(end)=='\' || filepath(end)=='/',
    filepath=filepath(1:end-1);
end;
%--------------------------------------------------------------------------

filename = dir([filepath, '/', old_base, '*', file_ext]);
filename = {filename.name};
[fnum,ind] = sort(str2double(cellfun(@(x) x(length(old_base)+1 : end-length(file_ext)), ...
    filename, 'UniformOutput', false)));
filename = filename(ind);
nn = length(fnum);
if nn==0,
    error('No specified files found in the directory!');
end;

% new file number
if keep_spacing,
    fnum = fnum-fnum(1)+numfmt(1);
else
    fnum = (0:nn-1)+numfmt(1);
end;
ndigit = floor(log10(fnum(end)))+1;
if ndigit<numfmt(2),
    ndigit = numfmt(2);
end;
ndigit = num2str(ndigit);

if ispc,
    for kk=1:nn,
        newfilename=[new_base, sprintf(['%0' ndigit 'd'],fnum(kk)), file_ext];
        dos(['rename "' filepath, '\', filename{kk}, '" "', newfilename, '"']);
    end;
else
    for kk=1:nn,
        newfilename=[new_base, sprintf(['%0' ndigit 'd'],fnum(kk)), file_ext];
        movefile([filepath,'/',filename{kk}],[filepath,'/',newfilename]);
    end;
end;
