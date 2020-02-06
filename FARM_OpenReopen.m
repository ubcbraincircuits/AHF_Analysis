function expe = FARM_OpenReopen(varOpen)

disp(['--- RELOAD DATA FROM: ' varOpen.filesave ': ---'])
for i = 1:10000
    filesave_i = [varOpen.filesave '_' num2str(i) '.mat'];
    try, 
        load(filesave_i); disp([' > load ' filesave_i]); expe(i) = expei; clear expei
    end
end, 
    