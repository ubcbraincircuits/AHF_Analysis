% function FARM_seed(varOpen)
% Function who seed the calcium activity and compare it to behavior
%    FARM_Preprocess : Function who load the data on the workspace
%    FARM_Process : Function who process the data before exploration
%  * FARM_seed : Function who seed the calcium activity and compare it to behavior
%    FARM_video : Function who play or record video of dual calcium and behavior
%    FARM_averaging : Function who average trial and seed
%
% INPUT: varOpen.Seed. (Structure)
%    .n: (scalar) number of behavior roi (default: 2)
%    .redo_roi: (scalar, 0/1) redraw the roi of behavior (default: 0)
%    .blue_correction: (scalar 0/1) display the reflectance (default: 1)
%    .smooth_beh: (scalar, in sec) smooth the temporal profile of behavior (default: .5)
%    .norm_beh: (scalar, 0/1) normalize behavior quantification from 0 to 1 (default: 0)
%    .norm_HyCM: (scalar, 0/1) normalize behavior quantification by the average (default: 0)
%    .sameRange_HyCM: (scalar, 0/1) normalize hybrid correlation map by average (default: 0)
%    .roimap: (scalar) map used for mapping (1: RMS/SD mixt, 2: [TO DEFINE] default: 1)
%    .win: (2 scalars [min max] or [] for all) window of display (default: [])
%    .jitter: (scalar, in sec) lag jitter of calcium relative to behavior (default: 0)
%    .coef: (scalar) exponent of the image behavior display IMG^coef (default: .5)
%    .ratio_ref: (scalar, 0 to 1) proportion of background image in behavior display (default: .5)
%    .k: (scalar) size of the average roi of clicking (square of (k*2)+1 edge) 0=1x1, 1=3x3, 2=5x5 (default: 1)
%    .ShowTask: (scalar, 0/1) display the task text (default: 0)
%    .showlist: (array of cell) boundaries of Allen Inst. areas: A AC AL AM AU BC FL HL L LI M1 M2 MO NO PL PM POR RL RS S2 TEA TR UN V1
%                (default: {'V1','HL','FL','BC','M1','RS'})
%
%         NB: to remove any field of a structure, use rmfield:
%             e.g. varOpen = rmfield(varOpen,'Seed'); to use only default
%             e.g. varOpenv.Seed= rmfield(varOpenv.Seed,'n'); to only remove one sub-field

function [varOpen,expe] = FARM_seed_nonc(varOpen,expe)



% -------------------- 1. Define default variable -------------------------
disp('---------------------------------------------------------------------')
try, disp(['              varOpen.Seed.blue_correction=' num2str(varOpen.Seed.blue_correction)]); 
    catch, varOpen.Seed.blue_correction = 1; disp('(not defined) varOpen.Seed.blue_correction=1'); end % set the default reflectance correction on
try, disp(['              varOpen.Seed.smooth_beh=' num2str(varOpen.Seed.smooth_beh)]); 
    catch, varOpen.Seed.smooth_beh = .5; disp('(not defined) varOpen.Seed.smooth_beh=.5'); end  % set the default behavior video binning at 0.5x
try, disp(['              varOpen.Seed.n=' num2str(varOpen.Seed.n)]); 
    catch, varOpen.Seed.n = 2; disp('(not defined) varOpen.Seed.n=2'); end  % set the default number of bhavior roi at 2
try, disp(['              varOpen.Seed.norm_beh=' num2str(varOpen.Seed.norm_beh)]); 
    catch, varOpen.Seed.norm_beh = 0; disp('(not defined) varOpen.Seed.norm_beh=0'); end  % set the default normalization of behavior to 0 (keep original range of absolute gradient)
try, disp(['              varOpen.Seed.norm_HyCM=' num2str(varOpen.Seed.norm_HyCM)]); 
    catch, varOpen.Seed.norm_HyCM = 0; disp('(not defined) varOpen.Seed.norm_HyCM=0'); end  % set the default normalization of hybrid correlation to 0 (not normalized by the average of every hybrid correlation map)
try, disp(['              varOpen.Seed.sameRange_HyCM=' num2str(varOpen.Seed.sameRange_HyCM)]); 
    catch, varOpen.Seed.sameRange_HyCM = 0; disp('(not defined) varOpen.Seed.sameRange_HyCM=0'); end % set the default value for normalizing the hybrid correlation at the same range    
try, disp(['              varOpen.Seed.redo_roi=' num2str(varOpen.Seed.redo_roi)]); 
    catch, varOpen.Seed.redo_roi = 0; disp('(not defined) varOpen.Seed.redo_roi=0'); end % set the default value to redo the behavior roi at each run (=0, re-use it when drawn)
try, disp(['              varOpen.Seed.roimap=' num2str(varOpen.Seed.roimap)]); 
    catch, varOpen.Seed.roimap = 1; disp('(not defined) varOpen.Seed.roimap=1'); end   % set the default value to define what behavior reference image to use (1=gradient+RMS) 
try, disp(['              varOpen.Seed.win=' num2str(varOpen.Seed.win)]); 
    catch, varOpen.Seed.win = []; disp('(not defined) varOpen.Seed.win=[]'); end % set the default time window (= full value)
try, disp(['              varOpen.Seed.jitter=' num2str(varOpen.Seed.jitter)]); 
    catch, varOpen.Seed.jitter = 0; disp('(not defined) varOpen.Seed.jitter=0'); end  % set the default jitter time use (=0)     
try, disp(); 
        for i = 1:length(varOpen.Seed.showlist), disp(['              varOpen.Seed.showlist=(' num2str(i) ')' char(varOpen.Seed.showlist(i))]); end % set the roi atlas border to show
    catch, varOpen.Seed.showlist = {'V1','HL','FL','BC','M1','RS'}; disp('(not defined) varOpen.Seed.showlist={''V1'',''HL'',''FL'',''BC'',''M1'',''RS''}'); end     
try, disp(['              varOpen.Seed.coef=' num2str(varOpen.Seed.coef)]); 
    catch, varOpen.Seed.coef = .5; disp('(not defined) varOpen.Seed.coef=.5'); end % set the default exponent for display, = image^0.5 = sqrt(image) to let more contribution of low values
try, disp(['              varOpen.Seed.ratio_ref=' num2str(varOpen.Seed.ratio_ref)]); 
    catch, varOpen.Seed.ratio_ref = .5; disp('(not defined) varOpen.Seed.ratio_ref=.5'); end % set the default proportion of reference image (50%)            
try, disp(['              varOpen.Seed.k=' num2str(varOpen.Seed.k)]); 
    catch, varOpen.Seed.k = 1; disp('(not defined) varOpen.Seed.k=1'); end  % set the default size of seed pixel, =1 (3x3)          
try, disp(['              varOpen.Seed.ShowTask=' num2str(varOpen.Seed.ShowTask)]); 
    catch, varOpen.Seed.ShowTask = 0; disp('(not defined) varOpen.Seed.ShowTask=0'); end   % don't show the text file for each task
try, disp(['              varOpen.Seed.colour_bar=' num2str(varOpen.Seed.colour_bar)]);
    catch, varOpen.Seed.colour_bar = []; disp('(not defined) varOpen.Seed.colour_bar=[]'); end 

%if varOpen.CBV == 0, varOpen.Seed.blue_correction = 0; disp('no reflectance channel, so varOpen.Seed.blue_correction=0'); end;
    
run_list = varOpen.varProc.run_list; 
range = [];

Tmin = min(varOpen.Seed.win); if length(Tmin)==0, Tmin=0; end; % offset starting time, =zero if no window
if varOpen.varOutput.sBe == 0 % if no behavior used
    varOpen.Seed.roimap = 0;% then, no roimap use to define any roi
end

if length(varOpen.Seed.win)==0, 
    varOpen.Seed.win = [0 varOpen.varProc.events(max(run_list)).time_window(2)];  % if no window define, use the full sequence
        use_expe = run_list ;
        
else 
    for i = 1:length(varOpen.Seed.win) % find the range of where the expes lie 
        for j = 1:length(run_list) % take min and max of seed window and see if it lies in expe by looking at time window for each expe 
            range_find = find((varOpen.varProc.events(run_list(j)).time_window(1) <= varOpen.Seed.win(i)) & (varOpen.varProc.events(run_list(j)).time_window(2)>= varOpen.Seed.win(i))); 
            if length(range_find) ~= 0 % from range defined in process, find the expes where these times lie 
                range = [range run_list(j)]; 
                break 
            else 
            end 
        end 
    end 
    use_expe = []; % run each expe in the specified min:max range and create array showing the expes are are interested in 
    for i = range(1):range(2) 
        expe_find = findstr(i,run_list); % find out which expes are part of run list and spit out array to take into account the following expes 
        if length(expe_find) ~= 0 
            use_expe = [use_expe i]; 
        else 
        end 
    end 
end 
    b_SF =[];
    b_SF_c= [];
    bxx = []; % time for behavioural data
    cxx = []; % calcium data 
    cat_beh = []; % concatenated behavioural frames 
    cat_cal = []; % concatenated calcium frames 
    for i = 1:length(use_expe)
        frame_n = []; % behavioural data for a single expe 
        frame_cal = []; % calcium data for a single expe 
        cx= []; % x time points for calcium data 
        bx = []; % x time corresponding behavioural  
        if i == 1 % for the first expe run, we take into account all the frames after the frame corresponding to the min window 
            frame_n = expe(use_expe(i)).Beh_cam(:,:,round((varOpen.Seed.win(1)-varOpen.varProc.events(use_expe(i)).time_window(1))*expe(use_expe(i)).realSF_cam/varOpen.varProc.tb+1): size(expe(use_expe(i)).Beh_cam,3)); 
            bx = linspace(varOpen.Seed.win(1),varOpen.varProc.events(use_expe(i)).time_window(2),size(frame_n,3)); 
            frame_cal = expe(use_expe(i)).If(:,:,round((varOpen.Seed.win(1)-varOpen.varProc.events(use_expe(i)).time_window(1))*expe(use_expe(i)).realSF/varOpen.varProc.tb+1):size(expe(use_expe(i)).If,3),:); 
            cx = linspace(varOpen.Seed.win(1),varOpen.varProc.events(use_expe(i)).time_window(2),size(frame_cal,3)); 
        elseif i == length(use_expe) % last expe run, where we take into account all frames that come before the max window 
            frame_n = expe(use_expe(i)).Beh_cam(:,:,1:round((varOpen.Seed.win(2)-varOpen.varProc.events(use_expe(i)).time_window(1))*expe(use_expe(i)).realSF_cam/varOpen.varProc.tb)); 
            bx = linspace(varOpen.varProc.events(use_expe(i)).time_window(1),varOpen.Seed.win(2),size(frame_n,3)); 
            frame_cal = expe(use_expe(i)).If(:,:,1:round((varOpen.Seed.win(2)-varOpen.varProc.events(use_expe(i)).time_window(1))*expe(use_expe(i)).realSF/varOpen.varProc.tb),:);
            cx = linspace(varOpen.varProc.events(use_expe(i)).time_window(1),varOpen.Seed.win(2),size(frame_cal,3)); 
        else % if it falls in between the first and last entries in the use_expe array, then take all frames into account 
            frame_n = expe(use_expe(i)).Beh_cam(:,:,:) ; 
            bx = linspace(varOpen.varProc.events(use_expe(i)).time_window(1), varOpen.varProc.events(use_expe(i)).time_window(2), size(frame_n,3)); 
            frame_cal = expe(use_expe(i)).If(:,:,:,:); 
            cx = linspace(varOpen.varProc.events(use_expe(i)).time_window(1), varOpen.varProc.events(use_expe(i)).time_window(2),size(frame_cal,3)); 
        end 
        
        cat_beh = cat(3,cat_beh,frame_n) ; % concatenate the behavioural data for all interested expes 
        bxx = [bxx bx] ; % x values corresponding to each frame (calculated via its individual real sampling frequencies) 
        varOpen.varProc.events(use_expe(i)).cx = cx ;
        varOpen.varProc.events(use_expe(i)).bx = bx ; % store individual time stamps of behavioural data
        
        cat_cal = cat(3,frame_cal,cat_cal) ; % concatenate the calcium data for each expe 
        cxx = [cxx cx]; % calculated array for the x values corresponding to each frame with respect to its specified sampling frequency 
    end

Beh_temp_file = [varOpen.working_folder 'Beh_temp.mat'];
try
    load(Beh_temp_file)  % if roiBeh defined, reload.
catch
    varOpen.Seed.redo_roi = 1; % instead, tag to redoit
end


couleur = imresize(jet,[varOpen.Seed.n 3],'nearest'); % LUT of color for behavior

q = who('JJJ'); % check if there is a map of K-means (behavior)
q2 = who('HyCMinv_RGB'); % check if there is a map of iinverted hybrid cross-correlation
q = length(q)+length(q2); % if none of the K-means of inverted hybrid cross-correlation is present, then map using map of movement
if q==0, roimap=1; end

if length(varOpen.Seed.showlist)>0 % if roi borders shown
    load Ai % load Allen Institute Data from toolbox
    clear listeroi, for i = 1:length(atlas.areatag) % for each area of Ai define
        listeroi(i) = 0; for j = 1:length(varOpen.Seed.showlist) % for each roi selected
            if strcmp(char(atlas.areatag(i)),char(varOpen.Seed.showlist(j)))==1, listeroi(i) = 1; end, end, end % if present, list if to show
    CA = atlas.CA(:,:,listeroi==1); % only keep the coordinates of selected borders
    BX = atlas.BX; % report bregma accordingly to bregma
    BY = atlas.BY+200;
    CA2 = zeros(size(CA,1)+200,size(CA,2),size(CA,3));
    CA2(201:200+size(CA,1),:,:)=CA;
    xx = round([BX-(varOpen.varOutput.L/2)*atlas.mmpp BX+(varOpen.varOutput.L/2)*atlas.mmpp]); 
    yy = round([BY-(varOpen.varOutput.L/2)*atlas.mmpp BY+(varOpen.varOutput.L/2)*atlas.mmpp]); 
    CA2 = CA2(yy(1):yy(2),xx(1):xx(2),:);
    B = [];
    for i = 1:size(CA2,3)
       BWCA = bwlabel(CA2(:,:,i));
       for j = 1:max(BWCA(:))
           B = [B bwboundaries(BWCA==j)];
       end
    end
end

if varOpen.Seed.roimap == 1
    RGB_REF = OIA_color(varOpen.varOutput.MM,gray);
    SDn = OIA_n(varOpen.varOutput.SD./varOpen.varOutput.MM).^varOpen.Seed.coef;
    Gn = OIA_n(varOpen.varOutput.G./varOpen.varOutput.MM).^varOpen.Seed.coef;
    RMSn = OIA_n(varOpen.varOutput.RMS./varOpen.varOutput.MM).^varOpen.Seed.coef;
    
    clear RGB
    RGB(:,:,1) = 255*OIA_n(SDn);
    RGB(:,:,3) = 255*OIA_n(Gn);
    RGB(:,:,2) = 255*OIA_n(RMSn);    
    
    RGB = uint8(RGB);
    RGB = RGB*(1-varOpen.Seed.ratio_ref) + RGB_REF*varOpen.Seed.ratio_ref;
    try, 
        RGB = insertText(RGB,[0 0],'SD','Font', 'Arial','FontSize',8,'TextColor','red','BoxColor','black','BoxOpacity',0);
        RGB = insertText(RGB,[0 8],'Grad','Font', 'Arial','FontSize',8,'TextColor','green','BoxColor','black','BoxOpacity',0);
        RGB = insertText(RGB,[0 16],'RMS','Font', 'Arial','FontSize',8,'TextColor','blue','BoxColor','black','BoxOpacity',0);
    end
elseif varOpen.Seed.roimap == 2
    RGB = JJJ;
    try
        RGB = insertText(RGB,[0 0],'behav K-means','Font', 'Arial','FontSize',8,'TextColor','white','BoxColor','black','BoxOpacity',0);
    end
elseif varOpen.Seed.roimap == 3
    RGB = HyCMinv_RGB;
    try,
        RGB = insertText(RGB,[0 0],'Hybrid Correlation','Font', 'Arial','FontSize',8,'TextColor','white','BoxColor','black','BoxOpacity',0);
    end
end

if varOpen.Seed.roimap ~= 0
    if varOpen.Seed.redo_roi == 1 % if define roi
        subplot(111)
        clear roiBeh % erase behavior roi
        for i = 1:varOpen.Seed.n % for each roi
            roiBeh(:,:,i) = OIA_roidraw(RGB); % draw roi
            roiBeh_disp = imdilate(roiBeh(:,:,i),ones(3)) - roiBeh(:,:,i); % define contour for display
            for j = 1:3
                RGB(:,:,j) = RGB(:,:,j).*uint8(~roiBeh_disp) + uint8(roiBeh_disp*256*couleur(i,j));
            end
        end
        save(Beh_temp_file,'roiBeh')
    else
        for i = 1:varOpen.Seed.n
            roiBeh_disp = imdilate(roiBeh(:,:,i),ones(3)) - roiBeh(:,:,i);
            for j = 1:3
                RGB(:,:,j) = RGB(:,:,j).*uint8(~roiBeh_disp) + uint8(roiBeh_disp*256*couleur(i,j));
            end
        end
    end
end


% 3. process Hybrid Correlation behavior to calcium
if varOpen.Seed.roimap ~= 0
    disp('Processing Hybrid Correlation ')
    JJ1d = single(reshape(cat_beh, size(cat_beh,1)*size(cat_beh,2), size(cat_beh,3))); % already extracted trials that we wanted 
    for i = 1:length(run_list) 
        Beh_cam = expe(run_list(i)).Beh_cam; % take individual behavioural data and reshape it so we can manipulate time dim later for correlation
        varOpen.Seed_expe(run_list(i)).JJ1ds = single(reshape(Beh_cam, size(Beh_cam,1)*size(Beh_cam,2), size(Beh_cam,3))); 
    end 
    clear be mini_be maxi_be% vector of mvt (roi x time)
    
    for i = 1:varOpen.Seed.n % for each roi
        roiBeh_i = roiBeh(:,:,i); % cut the roi
        for j = 2:size(JJ1d,2), % process of the density of mvt (absolute gradient) for each seed (to be used in averaging) 
            be(i,j) = mean(abs(single(JJ1d(roiBeh_i(:)==1,j)-JJ1d(roiBeh_i(:)==1,j-1))));
        end;
        if varOpen.Seed.smooth_beh>0, 
            be(i,:) = smooth(be(i,:),(varOpen.Seed.smooth_beh*(varOpen.SF/varOpen.varProc.tb))); 
        end
        %be(i,:) = OIA_n(be(i,:)); 
        mini_be(i)=min(be(i,:)); % min and max of each of the chosen ROIs
        maxi_be(i)=max(be(i,:));
    end
    
    % scale to all the roi
    mini_be = min(mini_be); % find the max and min between the 3 ROIs chosen 
    maxi_be = max(maxi_be);
 
    for i = 1:varOpen.Seed.n
        be(i,:) = be(i,:) - mini_be; 
        be(i,:) = be(i,:) / (maxi_be-mini_be); % spits out gradient that uses the range of max and min values 
    end  
    
    if varOpen.Seed.norm_beh == 1 % if normalized behavior
        be_average = mean(be,1); % average behavior activity between each roi
        for i = 1:varOpen.Seed.n
            be(i,:) = be(i,:) - be_average; % substraction of this average activity
        end;
    end

    jitter = round(varOpen.Seed.jitter*(varOpen.SF/varOpen.varProc.tb)); % vector jitter is local
    if jitter ~= 0 % if jitter present
        be = circshift(be,[0 jitter]); % behavior jittered befor hybrid corr
    end;
    for r = 1:length(run_list) % process each movement gradients individually for concatenation with correct SF 
        clear be_expe 
        for i = 1:varOpen.Seed.n 
            roiBeh_i = roiBeh(:,:,i); % each ROi 
            for j = 2:size(varOpen.Seed_expe(run_list(r)).JJ1ds,2) 
                be_expe(i,j) = mean(abs(single(varOpen.Seed_expe(run_list(r)).JJ1ds(roiBeh_i(:) ==1,j)-varOpen.Seed_expe(run_list(r)).JJ1ds(roiBeh_i(:)==1,j-1)))); 
            end 
            if varOpen.Seed.smooth_beh > 0
                be_expe(i,:)= smooth(be_expe(i,:),(varOpen.Seed.smooth_beh*(varOpen.SF/varOpen.varProc.tb))); 
            end
        mini_be_expe(i)=min(be_expe(i,:)); 
        maxi_be_expe(i) = max(be_expe(i,:)); 
        end 
        mini_be_expe = min(mini_be_expe); 
        maxi_be_expe = max(maxi_be_expe); 
        
        for i = 1:varOpen.Seed.n 
            be_expe(i,:) = be_expe(i,:) - mini_be_expe; 
            be_expe(i,:) = be_expe(i,:) / (maxi_be_expe-mini_be_expe) ; 
        end 
        
        if varOpen.Seed.norm_beh == 1 
            be_average_expe = mean(be_expe,1); 
            for i = 1:varOpen.Seed.n 
                be_expe(i,:) = be_expe(i,:) - be_average_expe; 
            end 
        end 
        if jitter ~= 0 
            be_expe = circshift(be_expe,[0 jitter]); 
        end 
        varOpen.Seed_expe(run_list(r)).be_expe = be_expe ; % save individual behavioural data in field (need for averaging) 
    end 

    if size(cat_cal,4)==1, warning('BLUE CORRECTION IMPOSSIBLE, NO CHANNEL PRESENT'); varOpen.Seed.blue_correction = 0; end
    HyCM = zeros(size(cat_cal,1),size(cat_cal,2),varOpen.Seed.n); % calculate hybrid correlation
    
    refl_fluo = []; 
    for i = 1:size(cat_cal,4)
        refl_fluo_i= []; 
        for j = 1:numel(use_expe) % Take calcium and behavioural data from each expe to be interpolated to have the same time points (same frames) and then interpolated so it can be correlated.
            if j == 1 
                I = expe(use_expe(j)).If(:,:,round((varOpen.Seed.win(1)-varOpen.varProc.events(use_expe(j)).time_window(1))*expe(use_expe(j)).realSF/varOpen.varProc.tb+1):size(expe(use_expe(j)).If,3),i); 
                Be = expe(use_expe(j)).Beh_cam(:,:,round((varOpen.Seed.win(1)-varOpen.varProc.events(use_expe(j)).time_window(1))*expe(use_expe(j)).realSF_cam/varOpen.varProc.tb+1):size(expe(use_expe(j)).Beh_cam,3)); 
            elseif j == numel(use_expe) 
                I = expe(use_expe(j)).If(:,:,1:(round((varOpen.Seed.win(2)-varOpen.varProc.events(use_expe(j)).time_window(1))*expe(use_expe(j)).realSF/varOpen.varProc.tb)),i); 
                Be = expe(use_expe(j)).Beh_cam(:,:,1:round((varOpen.Seed.win(2)-varOpen.varProc.events(use_expe(j)).time_window(1))*expe(use_expe(j)).realSF_cam/varOpen.varProc.tb)); 
            else 
                I = expe(use_expe(j)).If(:,:,:,i); 
                Be = expe(use_expe(j)).Beh_cam(:,:,:); 
            end 
            cx = varOpen.varProc.events(use_expe(j)).cx;
            bx = varOpen.varProc.events(use_expe(j)).bx; 

%     interpolate data trial by trial and concatenate; calcium will have
%     the same time points as the behavioural data by the end 
    
                [Xb,Yb,Zb] = meshgrid(1:size(I,1), 1:size(I,2),1:size(Be,3)) ; 
                [Xc,Yc,Zc] = meshgrid(1:size(I,1), 1:size(I,2),1:size(I,3)); 
                
                for c1 = 1:size(I,1) 
                    for c2 = 1:size(I,2) 
                        Zc(c1,c2,:) = cx; % assign time points to every frame
                        Zb(c1,c2,:) = bx; 
                    end 
                end 
                
                I= interp3(Xc,Yc,Zc,I,Xb,Yb,Zb); % interpolate it so I will have time points Zb 
                refl_fluo_i = cat(3,refl_fluo_i,I); % concatenated calcium frames with behaviour time points in time dim 
        end 
        refl_fluo = cat(4,refl_fluo,refl_fluo_i); % concatenate the refl and fluo data 
    end 
                
    for i = 1:varOpen.Seed.n % for each roi
        for ii = 1:size(refl_fluo,1) % for each brain pixel
            for jj = 1:size(refl_fluo,2)
                  s = single(squeeze(refl_fluo(ii,jj,1:size(refl_fluo,3),1))); % take all the frames corresponding to fluo of specified pixel 
                s(isnan(s)) = 0; % take all NaNs created from interpolation and replace with 0 to do correlation
                R=corrcoef(s,be(i,:)); HyCM(ii,jj,i)=R(1,2); % correlation coefficient; spits out matrix  
            end
        end
        HyCM(:,:,i) = HyCM(:,:,i).*imresize(varOpen.varOutput.roinan,size(HyCM(:,:,i))); % apply the mask
    end

    if varOpen.Seed.sameRange_HyCM == 1 % if normalization of hybrid corr maps from 0 to 1
        for i = 1:size(HyCM,3)
            HyCM(:,:,i) = HyCM(:,:,i) - min(min(HyCM(:,:,i)));
            HyCM(:,:,i) = HyCM(:,:,i) ./ max(max(HyCM(:,:,i)));
        end
    end

    if varOpen.Seed.norm_HyCM == 1 % if normalization of hybrid corr maps relatively to the average
        HyCM_average = mean(HyCM,3);
        for i = 1:varOpen.Seed.n
           HyCM(:,:,i) = HyCM(:,:,i) - HyCM_average;
        end;
    end

    HyCMd = []; % for display
    for i = 1:varOpen.Seed.n
        HyCMd = cat(2,HyCMd,HyCM(:,:,i));
    end;
    
    JJsb = imresize(single(cat_beh),1/4); % behavior 4x4 binned
    Gsb = zeros(size(JJsb)); % calculate the gradient 4x4 spatially binned
    for i = 2:size(JJsb,3)
        Gsb(:,:,i) = abs(JJsb(:,:,i)-JJsb(:,:,i-1));
    end
    
end

%% 4. clicking loop to associate calcium and behavior

% use concatenated behavioural and calcium data for the plot 

LineWidthValue = 1.5;
k = varOpen.Seed.k; % size of roi
q=who('x3'); if length(q)==0, x3=k+1;y3=k+1; button = 1 ; end % automatically set first position if not defined; x3, y3 and button defined by myginput
while button ~= 3 
    try, % automatically reset first position if not defined within roi
        ca = squeeze(mean(mean(refl_fluo(y3-k:y3+k,x3-k:x3+k,:,:),2),1)); % fluo data of all frames from pixel chosen with size of roi taken into account
    catch
        x3=k+1;y3=k+1;     
    end

    ca = squeeze(mean(mean(refl_fluo(y3-k:y3+k,x3-k:x3+k,:,:),2),1)); % fluo   
    
    if varOpen.Seed.blue_correction ~= 0, 
        re = ca(:,2); % refl
        ca = ca(:,1);
    end; 
    
%     bx = Tmin+[1:length(ca)]'/(varOpen.SF/varOpen.varProc.tb);
 
    fig1=figure(1);
    fig1.Renderer='Painters';
    subplot(212), 

    
    if varOpen.Seed.blue_correction ~= 0, 
        plot(bxx,ca-varOpen.Seed.blue_correction.*re,'k','LineWidth',LineWidthValue); hold on; % subtraction
        plot(bxx,varOpen.Seed.blue_correction.*re,'b','LineWidth',LineWidthValue); % reflectance 
        plot(bxx,ca,'g','LineWidth',LineWidthValue); % fluorescence 
    else
        plot(bxx,ca,'g','LineWidth',LineWidthValue); hold on; 
    end
    
    range_calcium = max(ca)-min(ca);
    base_beh = min(ca)+range_calcium;
    
    if varOpen.Seed.roimap ~= 0        
        for i = 1:varOpen.Seed.n
            ben = OIA_n(be(i,:)) * range_calcium;
            plot(bxx,base_beh+ben,'Color',couleur(i,:),'LineWidth',LineWidthValue); % behavioural data for every ROI
        end
        base_beh = base_beh + max(ben(:));
    end
    
    for i=1:length(run_list) % shows the LEDOFF times for each expe and indicates potential overlap
        plot([varOpen.varProc.events(run_list(i)).time_window(2) varOpen.varProc.events(run_list(i)).time_window(2)],[min(ca) base_beh],'r')
    end 
    
    plot(varOpen.varOutput.event_list,ones(size(varOpen.varOutput.event_list))*base_beh,'r.'); % licks 
    clear LR LRx
    for i = 2:length(varOpen.varOutput.event_list) % instantaneous licks and its curve fitting line 
        LR(i-1) = 1/(varOpen.varOutput.event_list(i)-varOpen.varOutput.event_list(i-1));
        LRx(i-1) = varOpen.varOutput.event_list(i); 
    end
    
    try,    
        plot(LRx,base_beh+LR,'c.');
        LRi = interp1(LRx,LR,bxx); % interpolation of the instantaneous licks calculated above
        LRi = smooth(LRi,(varOpen.Seed.smooth_beh*(varOpen.varOutput.SFr)));    
        hold on, plot(bxx,base_beh+LRi,'m-'); % curve fitting line for the instaneous licks 
        add_text = '/ corr slick vs. (roi): ';
        if varOpen.varOutput.sBe ~= 0
            for i = 1:varOpen.Seed.n 
                RR = corrcoef(LRi,be(i,:)); % correlation between instaneous licks and behavioural gradient 
                add_text = [add_text ['(' num2str(i) ')r=' num2str(myround(RR(1,2),2)) ' ']];
            end;
        end
    end
    for i= 1:length(varOpen.varOutput.event_list_task_info)
        if (length(strfind(char(varOpen.varOutput.event_list_task_info(i)),'GO=2'))|length(strfind(char(varOpen.varOutput.event_list_task_info(i)),'GO=-2'))|length(strfind(char(varOpen.varOutput.event_list_task_info(i)),'GO=-4'))) ==1
            plot(varOpen.varOutput.event_list_task(i),base_beh,'gv'); % indicate the times of the go trials 
        else 
            plot(varOpen.varOutput.event_list_task(i), base_beh, 'rv'); % indicate no go trials 
        end
    end
    
    plot(varOpen.varOutput.event_list_reward,ones(size(varOpen.varOutput.event_list_reward))*base_beh,'b*'); % reward 
    if varOpen.Seed.ShowTask == 1 % show the details of each GO trial where all triangles are marked 
        for i = 1:length(varOpen.varOutput.event_list_task)
            text(double(varOpen.varOutput.event_list_task(i)),double(base_beh),char(varOpen.varOutput.event_list_task_info(i)),'Rotation',45,'Color',[0 1 0])
        end
    end
    
    
    xlabel('time (s)'); ylabel('DF/F (%)'); hold off; % display axis and legend eventually
    if varOpen.Seed.blue_correction ~= 0, legend('subs',['refl x' num2str(varOpen.Seed.blue_correction)],'fluo'); end
    xlim([min(bxx) max(bxx)]); % window of time shown
    
    
    REF_disp = OIA_color(single(varOpen.varOutput.REF(:,:,1)),gray);

    CM = OIA_SA_seedpixel(refl_fluo(:,:,:,1),x3,y3);
    CM = CM.*imresize(varOpen.varOutput.roinan,size(CM));
    
    titre = 'Reflectance / Seed pixel Corr';
    disp_panel = 2;
    try,
        CM2 = OIA_SA_seedpixel(refl_fluo(:,:,:,2),x3,y3);
        CM2 = CM2.*imresize(varOpen.roinan,size(CM2));
        CM = cat(2,CM,CM2);
        disp_panel = disp_panel + 1;
        titre = 'Reflectance / Seed pixel Corr (Fluo/Refl)';
    end
    
    CM_disp = CM .*((CM>=0)&(CM<=1));
    CM_disp = OIA_color(CM,jet);
    Combined_disp = cat(2,REF_disp,CM_disp);
    
    if varOpen.Seed.roimap == 0, subplot(211); 
    else, subplot(232); end
    
    imshow(Combined_disp); colorbar; colormap jet
    title(titre)
    
    scaleBar = size(varOpen.varOutput.REF,1)/varOpen.varOutput.L;
    hold on, 
    for i = 1:disp_panel
        offset = size(CM,1)*(i-1);
        plot(x3+offset,y3,'ow'); 
        plot([size(varOpen.varOutput.REF,1)/20 scaleBar+size(varOpen.varOutput.REF,1)/20]+offset,[size(varOpen.varOutput.REF,1)/20 size(varOpen.varOutput.REF,1)/20],'w','LineWidth',2); 
        if length(varOpen.Seed.showlist)>0, % put this in function 
            for j = 1:length(B), 
                Bi = cell2mat(B(j))*size(CM,1)/length(BWCA); 
                offset_x = (size(CM,1)/2) - varOpen.varOutput.bx;
                offset_y = (size(CM,1)/2) - varOpen.varOutput.by;
                plot(Bi(:,2)-offset_y+offset,Bi(:,1)-offset_x,'w-'); 
            end; 
        end
        plot(varOpen.varOutput.by,varOpen.varOutput.bx,'xk'); 
        plot(varOpen.varOutput.by,varOpen.varOutput.bx,'+w'); 
    end
    hold off
    


        
    if varOpen.Seed.roimap ~= 0
        HyCMinv = zeros(size(JJsb,1),size(JJsb,2));
        for i = 1:size(JJsb,1)
            for j = 1:size(JJsb,2)
                bei = squeeze(squeeze(Gsb(i,j,:)));
                if varOpen.Seed.blue_correction ~= 0
                    val = corrcoef(bei,(ca-varOpen.Seed.blue_correction.*re));
                else
                    val = corrcoef(bei,ca);
                end
                HyCMinv(i,j) = val(1,2);
            end
        end
    

            
        range = myround([min(HyCMinv(:)) max(HyCMinv(:))],2);

        HyCMinv_RGB = OIA_color(HyCMinv,jet);
        HyCMinv_RGB = imresize(HyCMinv_RGB,4);

        LUT = imresize([range(1):(range(2)-range(1))/100:range(2)],[size(HyCMinv_RGB,1)/8 size(HyCMinv_RGB,2)/3]);
        LUT = OIA_color(LUT,jet);

        HyCMinv_RGB = HyCMinv_RGB*(1-varOpen.Seed.ratio_ref) + RGB_REF*varOpen.Seed.ratio_ref;
        HyCMinv_RGB(1:size(LUT,1),1:size(LUT,2),:)=LUT;
        try
            HyCMinv_RGB = insertText(HyCMinv_RGB,[0 0],['r=' num2str(range(1))],'Font', 'Arial','FontSize',12,'TextColor','white','BoxColor','black','BoxOpacity',0);
            HyCMinv_RGB = insertText(HyCMinv_RGB,[size(LUT,2)*.66 0],[num2str(range(2))],'Font', 'Arial','FontSize',12,'TextColor','white','BoxColor','black','BoxOpacity',0);
        end

        HyCMinv_RGB = cat(1,RGB,HyCMinv_RGB); 
        subplot(231)
        imshow(HyCMinv_RGB);
        title('Hybrid correlation (calcium to behav)')

 
        mini = nanmin(HyCMd(:))
        maxi = nanmax(HyCMd(:))
        superMax = max([abs(mini) abs(maxi)]);
        
        if length(varOpen.Seed.colour_bar)== 0 
        subplot(233); 
        imshow(HyCMd,[-superMax superMax]); 
        ligne = 1;
        else 
            subplot(233);
            imshow(HyCMd,[varOpen.Seed.colour_bar(1) varOpen.Seed.colour_bar(2)]);
            ligne = 1; 
        end 
            

        colorbar; colormap jet
        hold on, for i = 1:varOpen.Seed.n, plot(x3+size(refl_fluo,2)*(i-1),y3,'or'); end; hold off
        title('Hybride cross-corr (behav to calcium)')
        for i = 1:varOpen.Seed.n, 
            text((i-1)*size(cat_cal,1)+1,(ligne+0.1)*size(refl_fluo,1),['roi=' num2str(i)],'Color',couleur(i,:)); end

        for i = 1:varOpen.Seed.n
            if length(varOpen.Seed.showlist)>0, % module to display ABI atlas $$$$$$$$$$$$$$$$$$
                for j = 1:length(B), 
                    Bi = cell2mat(B(j))*size(HyCMd,1)/length(BWCA); 
                    offset_x = (size(HyCMd,1)/2) - varOpen.varOutput.bx;
                    offset_y = (size(HyCMd,1)/2) - varOpen.varOutput.by;
                    hold on, 
                    plot(Bi(:,2)-offset_y+(i-1)*size(HyCMd,1),Bi(:,1)-offset_x,'w-'); 
                    hold off, 
                end; 
            end 
        end
    end
    
    seed_temp_file = [varOpen.working_folder 'seed_temp.mat'];
    seed_info.ca = ca;
    seed_info.conc_cal = refl_fluo; % interpolated calcium data saved in temporary file 
    try, seed_info.re = re; end
    try, seed_info.be = be; end
    seed_info.bxx = bxx; % time points corresponding to behavioural and calcium data 
    seed_info.couleur = couleur;
    try, seed_info.roiBeh = roiBeh; end
    seed_info.x3 = x3;
    seed_info.y3 = y3;
    seed_info.base_beh = base_beh;
    save(seed_temp_file,'seed_info')
    
    [x3,y3,button]=myginput;
    while x3>size(refl_fluo,2)
        x3 = x3-size(refl_fluo,2);
    end
    while y3>size(refl_fluo,1)
        y3 = y3-size(refl_fluo,1);
    end
    
    
end
