function varOpen = FARM_OpenDefault(varOpen)

varOpen.file

disp('--------- SET DEFAULT VARIABLES --------')
try, disp(['              varOpen.mouse=' varOpen.mouse]); 
    catch, varOpen.mouse = ''; disp('(not defined) varOpen.mouse='''''); end

try, 
    disp(['              varOpen.file=' varOpen.file]); 
catch, 
    try
        for i = 1:length(varOpen.file)
            disp(['              varOpen.file(' num2str(i) ').rep=' varOpen.file(i).rep]); 
        end   
    catch
        varOpen.file = []; 
        disp('(not defined) varOpen.file=[]'); 
    end 
end
        
if length(varOpen.file)==0 disp('choose a file'); [file,path] = uigetfile('*.*'); varOpen.file=[path file]; end
    try, loca = findstr(varOpen.file(1).rep,'\'); varOpen.folder = varOpen.file(1).rep(1:loca(end));
    catch, loca = findstr(varOpen.file,'\'); varOpen.folder = varOpen.file(1:loca(end)); end  
try, disp(['              varOpen.filesave=' varOpen.filesave]); 
    catch, varOpen.filesave = [varOpen.folder 'dataset_' ]; disp(['(not defined) varOpen.filesave=' varOpen.filesave]); end
    varOpen.filesave = [varOpen.filesave varOpen.mouse]; 
try, disp(['              varOpen.SF=' num2str(varOpen.SF)]); 
    catch, varOpen.SF = 30; disp('(not defined) varOpen.SF=30'); end
try, disp(['              varOpen.BinningBeh=' num2str(varOpen.BinningBeh)]);  
    catch, varOpen.BinningBeh = .25; disp('(not defined) varOpen.BinningBeh=1/4'); end
try, disp(['              varOpen.reprocess=' num2str(varOpen.reprocess)]); 
    catch, varOpen.reprocess = 0; disp('(not defined) varOpen.reprocess=0'); end
try, disp(['              varOpen.format={''' char(varOpen.format) '''}']); 
    catch, varOpen.format =  {'raw','tif'}; disp('(not defined) varOpen.format={''raw'',''tif''};'); end
try, disp(['              varOpen.tb=' num2str(varOpen.tb)]); 
    catch, varOpen.tb = 1; disp('(not defined) varOpen.tb=1'); end
try, disp(['              varOpen.X=' num2str(varOpen.X)]); 
    catch, varOpen.X = 64; disp('(not defined) varOpen.X=64'); end
try, disp(['              varOpen.L=' num2str(varOpen.L)]);  
    catch, disp('(not defined) varOpen.L=10.25 if raw (or 8.6 if tif)'); end
try, disp(['              varOpen.CBV=' num2str(varOpen.CBV)]);  
    catch, disp('(not defined) varOpen.CBV will be requested'); end
    
try, disp(['              varOpen.registParam.f=' num2str(varOpen.registParam.f)]); 
    catch, disp('(not defined) varOpen.registParam.f=round([varOpen.X/20 varOpen.X/60])'); varOpen.registParam.f = round([varOpen.X/20 varOpen.X/60]); end
try, disp(['              varOpen.registParam.dx=' num2str(varOpen.registParam.dx)]);  
    catch, disp('(not defined) varOpen.registParam.dx=round(varOpen.X/12)'); varOpen.registParam.dx = round(varOpen.X/12); end
try, disp(['              varOpen.registParam.oo=' num2str(varOpen.registParam.oo)]);  
    catch, disp('(not defined) varOpen.registParam.oo=[-10:2:10]'); varOpen.registParam.oo=[-10:2:10]; end

try, disp(['              varOpen.filecam=' varOpen.filecam]); 
    catch, disp('(not defined) varOpen.filecam'); end    
try, disp(['              varOpen.filecam2=' varOpen.filecam2]); 
    catch, disp('(not defined) varOpen.filecam2'); end   
try, disp(['              varOpen.filecam3=' varOpen.filecam3]); 
    catch, disp('(not defined) varOpen.filecam3'); end   
 
% put a "\" at the end of folder in case if was forgotten
if isstruct(varOpen.file)==1, % if this is a list of file in structure
    if isdir(varOpen.file(1).rep)==1 % if this is a folder, I just check if there is a \ 
        for i = 1:length(varOpen.file)
            if varOpen.file(i).rep(end)~='\', 
                varOpen.file(i).rep = [varOpen.file(i).rep '\'];
            end
        end
    end
else, % if this is a single string file
    if isdir(varOpen.file)==1
        if varOpen.file(end)~='\', 
            varOpen.file = [varOpen.file '\'];
        end
    end
end

filesave_i = [varOpen.filesave '_1.mat'];
varOpen.test_filesave = dir(filesave_i);
if ((length(varOpen.test_filesave)~=0)&(varOpen.reprocess==1))
    disp(' -------------------------------------------------------------------')
    disp('  ')
    warning(['!!!!! you are going to overwrite your previous preprocessed data from: ' filesave_i '...'])
    q = input('     would you like to continue? (yes:1/no:0):');
    if q==0, varOpen.reprocess=0; end
end

if varOpen.reprocess == 1, 
    varOpen.test_filesave = []; 
    resu = dir([varOpen.filesave '*.mat']);
    rep = varOpen.filesave(1:max(findstr(varOpen.filesave,'\'))); 
    for i = 1:length(resu)
        delete([rep resu(i).name])
        disp([rep resu(i).name ' deleted !'])
    end
end

if isstruct(varOpen.file) == 0, % if only one file, convert into structure
    file2 = varOpen.file; 
    varOpen = rmfield(varOpen,'file');
    varOpen.file(1).rep = file2; 
end;

