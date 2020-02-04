% [expe,varOpen] = DualCaBe_Preprocess(varOpen); 
%    or [expe,varOpen] = DualCaBe_Preprocess; if manually defined
% Function who load the data on the workspace
%  * FARM_Preprocess : Function who load the data on the workspace
%    FARM_Process : Function who process the data before exploration
%    FARM_seed : Function who seed the calcium activity and compare it to behavior
%    FARM_video : Function who play or record video of dual calcium and behavior
%    FARM_averaging : Function who average trial and seed
%
% INPUT: varOpen (Structure)
%    .file: folder or file (or list) for calcium data
%            default: manual opening
%           EXAMPLE 1 (one file):
%               varOpen.file = 'D:\data_example_homecage_2017\M210608298_stim_1505656729.raw';
%           EXAMPLE 2 (n files):
%               varOpen.file(1).rep = 'D:\data_example_homecage_2017\M210608298_stim_1505656729.raw';
%               varOpen.file(2).rep = 'D:\data_example_homecage_2017\M210608298_stim_1505656992.raw';
%           EXAMPLE 3 (one folder): 
%               varOpen.file = 'D:\data_example_homecage_2017\';
%           EXAMPLE 4 (n files):
%               varOpen.file(1).rep = 'D:\data_example_homecage_2017\';
%               varOpen.file(2).rep = 'D:\data2_example_homecage_2018\';  
%           EXAMPLE 5 one file manually define):
%               varOpen.file = []; or not defined : for manually defined
%    .mouse: (string) mouse name
%            default: ''
%    .X: (scalar) size of the calcium data (after opening and spatial binning),
%            default: 64 (=64x64)
%    .tb: (scalar) temporal binning
%            default: 1 (= no binning)
%    .SF (scalar) sampling frequency, default: 30 
%    .format (cells) type of file used in folder, default : {'raw','tif'};
%    .L (scalar) width of the image, default 10.25 if raw, 8.6 if tif
%    .CBV (scalar 0 or 1) if CBV considerer (1=yes, 0=no), default: requested
%           if raw: keep green and blue frame
%           if tif: stronbing, interleave green reflectance and fluo
%    PARAMETER OF REGISTRATION IF MULTIPLE RUN:
%       .registParam.f (scalar), spatial bandpass to highlight vessels
%            default: round([varOpen.X/20 varOpen.X/60])
%       .registParam.dx (scalar), range of shift to test
%            default: round(varOpen.X/12)
%       .registParam.oo (scalar), range of angle to test
%            default: [-10:2:10]
%    .filecam:     behavior videos (view 1)
%       .filecam2: behavior videos (view 2)
%       .filecam3: behavior videos (view 3)
%             default: no behavior video (not defined)
%           EXAMPLE 1 (one file):
%               varOpen.filecam = 'D:\data_example_homecage_2017\210608298__1505656729.mp4';
%           EXAMPLE 2 (n files):
%               varOpen.filecam(1).rep = 'D:\data_example_homecage_2017\210608298__1505656729.mp4';
%               varOpen.filecam(2).rep = 'D:\data_example_homecage_2017\210608298__1505656992.mp4';
%           EXAMPLE 3 (one or multiple folder, match with calcium file name, in the same folder): 
%               varOpen.filecam = [];
%           EXAMPLE 4 (one single folder were all the videos are, independently of where the calcium file are): 
%               varOpen.filecam = 'D:\videotodelete\';
%    .BinningBeh (scalar), spatial binning of the behavior data, 
%             default = 1/4
%    .filesave (string) mat file where the data are stored for future re-opening
%             default: [folder]\dataset_[x].mat
%    .reprocess (scalar 0 or 1), reprocess of not .filesave mat file
%             default: 0
%
% OUTPUT: varOpen (updated)
%         expe Struture of data
%            .I: [X,Y,F,C] calcium data X,Y image for each frame F and
%                   channel C
%            .Is: average profile of I
%            .by: bregma x (pixel)
%            .bx: bregma y (pixel)
%            .L: width (mm)
%            .roi: roi mask
%            .SF: sampling frequency
%            .filename: original file name
%            .exten: extension 
%            .date: ASCII data
%            .datenum: real date
%            .d: index of the file if multiple folder
%            .be: structure of behavior data:
%                .be(i).videofilename: videofile name
%                .be(i).J: matrix of video [X,Y,F]
%                .be(i).s: average over time
%                .be(i).sJ_bin: binning of average over time
%                .be(i).sJ_after: period of calcium activity
%                .be(i).SD: standard deviation vs time
%                .be(i).MM: average vs time
%                .be(i).G: average of gradient vs time
%                .be(i).RMS: root mean square vs time

function [expe,varOpen] = DualCaBe_Preprocess(varOpen); 
    try, disp(['              varOpen.delete=' num2str(varOpen.delete)]);
        catch, varOpen.delete= 1; disp('(not defined) varOpen.delete=1'); end 
    
if nargin==0, varOpen.file = []; end

varOpen = FARM_OpenDefault(varOpen); % to set default variables if not defined

if length(varOpen.test_filesave)~=0 % if data already exist
    expe = FARM_OpenReopen(varOpen); % re-open if processed before   
else
    [expe,varOpen] = FARM_openraw(varOpen); % open calcium files (from raw or tif)
    expe = FARM_threshold(expe); % set the threshold of calcium files
    expe = FARM_opentextCage(expe); % associate each trial with task information from the text file
    expe = FARM_saturation(expe); % discard manually trial with too much saturated pixels or movements
    expe = FARM_registration(expe,varOpen); % automatic registration of calcium data
    expe = FARM_roi(expe,varOpen); % draw the roi and identify bregma

    
    try, disp(['Open video from ' varOpen.filecam]), % open behavior camera 1
        expe = FARM_behavVideo(expe,varOpen.filecam,varOpen,1); 
    end
    
    try, disp(['Open video from ' varOpen.filecam2]), % open behavior camera 2
        expe = FARM_behavVideo(expe,varOpen.filecam2,varOpen,2); 
    end
    
    try, disp(['Open video from ' varOpen.filecam3]), % open behavior camera 3
        expe = FARM_behavVideo(expe,varOpen.filecam3,varOpen,3); 
    end 
    
    FARM_save(expe,varOpen); % save data for future reopening
    
end   

FARM_ShowBehavior(expe); % show the SD and location of movement if behavior present
if varOpen.delete == 1
disp('Deleting previous temp files of seed analysis:')
delete([varOpen.folder 'Beh_temp.mat']) % delete temp file for behavior mapping
delete([varOpen.folder 'Ave_temp.mat']) % delete temp file for averaging mapping
else 
end 

