function [expe,varOpen] = FARM_openraw(varOpen)

try, varOpen.CBV; catch, 
    disp('----------------------------------------------------------------')
    disp('varOpen.CBV not defined')
    disp('    ?varOpen.CBV=1 reflectance recorded')
    disp('    ?varOpen.CBV=0 no reflectance')
    disp(' ')
    varOpen.CBV = input('Is there reflected light? (1/0)'); end


if isdir(varOpen.file(1).rep)==1 % if this is a folder
    disp('----------- [MULTIPLE FILES OPEN FROM ONE FOLDER (OR SERIES OF FOLDER)] -------------')
    ind = 0;
    for j = 1:length(varOpen.file) % for each folder listed

        % list all the file (in the current folder) with the extensions listed in varOpen.format
        resu = [];
        if length(varOpen.mouse) == 0, current_mouse = ''; else, current_mouse = [varOpen.mouse '*']; end
        for i = 1:length(varOpen.format) % for each extensions listed in varOpen.format
            disp(varOpen.format)
            resu = [resu; dir([varOpen.file(j).rep '*' current_mouse '.' char(varOpen.format(i))])]; % concatenate list of file corresponding      
        end
        for i = 1:length(resu) % for each file of each folder
            disp('Hello2')
            tic, ind = ind + 1;
            disp(['adding file (' num2str(i) '/' num2str(length(resu)) '): ' varOpen.file(j).rep resu(i).name])
            current_file = [varOpen.file(j).rep resu(i).name]; % name of the current file
            [I,L] = FARM_openraw_singleFile(current_file,varOpen); % open one single file
            % expe is the structure where all the data will be loaded
            try, expe(ind).L = varOpen.L; catch, expe(ind).L = L; end; % add the info about the width L
            expe(ind).I = I; % add the matrix
            expe(ind).d = j; % add the index of the folder
            resu2 = dir(current_file); % check the information about the file
            expe(ind).date = resu2.date; % report the date (ASCII)
            expe(ind).datenum = resu2.datenum; % report the date (real)
            expe(ind).filename = current_file; % report the name of the file
            disp(['added in t=' num2str(toc)]) % report the duration of this step
        end
    end
       
else % if this is a single file (or a serie of file)
    disp('----------- [SINGLE FILE OPEN (OR SERIES OF SINGLE FILE)] -------------')
    for j = 1:length(varOpen.file) % for each file listed
        disp(varOpen.file(j).rep)
        [I,L] = FARM_openraw_singleFile(varOpen.file(j).rep,varOpen);
        % expe is the structure where all the data will be loaded
        try, expe(j).L = varOpen.L; catch, expe(j).L = L; end
        expe(j).I = I; 
        expe(j).d = 0;
        resu2 = dir(varOpen.file(j).rep);
        expe(j).date = resu2.date;
        expe(j).datenum = resu2.datenum;
        expe(j).filename = varOpen.file(j).rep;
    end
end

