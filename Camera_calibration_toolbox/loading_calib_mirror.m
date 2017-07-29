% % Old_Multimirror_Calib_Results*.mat will not be loaded, you have to
% % load old calibration results manually if you want to use them.

dir('Multimirror_Calib_Results*.mat');

datasets = dir('Multimirror_Calib_Results*.mat');
Ndata = length(datasets);  % number of Multimirror_Calib_Results*.mat
if Ndata==0,
    fprintf(1,'Calibration results file not found in current directory!\n');
    return;
elseif Ndata>1,
    fprintf(1,['Multiple ''Multimirror_Calib_Results*.mat'' files found in current directory!\n' ...
        'You may merge these calibration results together if they are from the same camera!\n']);
    call_flag = input('Merge these calibration results into one or not? ([]=yes, other=no) ','s');
    call_flag = isempty(call_flag);
    if call_flag,
        merge_datasets_mirror;
        return;
    else
        fprintf(1,'You have to determine which calibration result file to load ...\n');
        n = 0;
        while ~n,
            n =  input(['Which file do you want to load? (choose from: 1~' num2str(Ndata) ', []=1)']);
            if isempty(n),
                n = 1;
            else
                n = round(n);
                if n<1 || n>Ndata,
                    fprintf(1,'Unexpected input!\n');
                    n = 0;
                end;
            end;
        end;
        save_name = datasets(n).name;
        fprintf(1,'\nLoading calibration results from %s ...\n',save_name);
    end;
else
    save_name = datasets(1).name;
    fprintf(1,'\nLoading calibration results from %s ...\n',save_name);
end;

load(save_name);
fprintf(1,'done\n');