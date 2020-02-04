% function FARM_average(varOpen)
% Function who average trial and seed
%    FARM_Preprocess : Function who load the data on the workspace
%    FARM_Process : Function who process the data before exploration
%    FARM_seed : Function who seed the calcium activity and compare it to behavior
%    FARM_video : Function who play or record video of dual calcium and behavior
%  * FARM_averaging : Function who average trial and seed
%
% INPUT: varOpen.ave. (Structure)
%    .qtype: (scalar) trigger mode (default: -1) = first roi
%             0: lick tag 
%             1: reward
%             2: GO=2
%             3: GO=-4
%  -[roi index]: profile of the roi (eg -2 for activity of roi #2)
%    .pre (scalar) pre-window in sec (default: 5)
%    .refr_period (scalar) refractory period in sec (only for manual trigger: qtype<0, default: 5)
%    .temp_filter (scalar) temporal fiter in sec (only for manual trigger (qtype<0) or mvt discard (mvt_roi>0) , default: 2)
%    .redo (scalar, 0/1) restart the process again (default: 0)
%    .range ([min max] or []) range displayed (default: [minimum maximum]
%    .mvt_roi = (scalar or 0) to discard trial if mvt within this roi is above threshold determined by mvt_kSD, see below (default = 0)
%    .mvt_kSD (real) to discard trial if mvt above mean+mvt_kSD*SD (default = 1)
%    .mvt_pre (real, 0 to 1) to calcule mvt on % of the window (default = 1)
%              e.g. .5: =50%, check mvt on the first half (-pre to 0), so preceeding trial tag only
%              e.g. 1 : =100%, check mvt on the full window from -pre to pre
%
%         NB: to remove any field of a structure, use rmfield:
%             e.g. varOpen = rmfield(varOpen,'ave'); to use only default
%             e.g. varOpenv.ave= rmfield(varOpenv.ave,'pre'); to only remove one sub-field

function [varOpen,expe] = FARM_averaging_teri_new(varOpen,expe)

disp('---------------------------------------------------------------------')
try, disp(['              varOpen.ave.qtype=' num2str(varOpen.ave.qtype)]); 
    catch, varOpen.ave.qtype = -1; disp('(not defined) varOpen.ave.qtype=-1'); end; % set the default value of tag use for averaging = 1 (=mvt evoked, roi 1)   
try, disp(['              varOpen.ave.pre=' num2str(varOpen.ave.pre)]); 
    catch, varOpen.ave.pre = 5; disp('(not defined) varOpen.ave.pre=5'); end; % set the default value of prestim baseline = 5s    
try, disp(['              varOpen.ave.refr_period=' num2str(varOpen.ave.refr_period)]); 
    catch, varOpen.ave.refr_period = 5; disp('(not defined) varOpen.ave.refr_period=5'); end; % set the default value of refractory period = 5s  
try, disp(['              varOpen.ave.temp_filter=' num2str(varOpen.ave.temp_filter)]); 
    catch, varOpen.ave.temp_filter = 2; disp('(not defined) varOpen.ave.temp_filter=2'); end; % set the default value of temporal filter to identify spikes = 2s 
try, disp(['              varOpen.ave.redo=' num2str(varOpen.ave.redo)]); 
    catch, varOpen.ave.redo = 0; disp('(not defined) varOpen.ave.redo=0'); end; % to redo the processing each time before clicking, = 0 (no)
try, disp(['              varOpen.ave.range=' num2str(varOpen.ave.range)]); 
    catch, varOpen.ave.range = []; disp('(not defined) varOpen.ave.range=[]'); end; % range of display (in DF/F %), default = [] = full range
try, disp(['              varOpen.ave.mvt_roi=' num2str(varOpen.ave.mvt_roi)]); 
    catch, varOpen.ave.mvt_roi = 0; disp('(not defined) varOpen.ave.mvt_roi=0'); end; % set the default value to discard trials that are above threshold in roi = 0
try, disp(['              varOpen.ave.mvt_kSD=' num2str(varOpen.ave.mvt_kSD)]); 
    catch, varOpen.ave.mvt_kSD = 1; disp('(not defined) varOpen.ave.mvt_kSD=1'); end; 
try, disp(['              varOpen.ave.mvt_pre=' num2str(varOpen.ave.mvt_pre)]); 
    catch, varOpen.ave.mvt_pre = 1; disp('(not defined) varOpen.ave.mvt_pre=1'); end;

Ave_temp_file = [varOpen.working_folder 'Ave_temp.mat'];

try
    load(Ave_temp_file)  % if trial contactenation file exit, reload.
    % this file contain 'I4' and original matrix: 'I2'
catch
    varOpen.ave.redo = 1; % If there's no averaging file, it will redo the whole process
end

seed_temp_file = [varOpen.working_folder 'seed_temp.mat']; % file of seed info (such as behavior)
load(seed_temp_file)
disp(['loaded from ' seed_temp_file])
%     seed_info.ca = ca;
%     seed_info.conc_cal = refl_fluo; % interpolated calcium data saved in temporary file 
%     try, seed_info.re = re; end
%     try, seed_info.be = be; end
%     seed_info.bxx = bxx; % time points corresponding to behavioural and calcium data 
%     seed_info.couleur = couleur;
%     try, seed_info.roiBeh = roiBeh; end
%     seed_info.x3 = x3;
%     seed_info.y3 = y3;
%     seed_info.base_beh = base_beh;
% 

run_list= varOpen.varProc.run_list ; 
if varOpen.ave.redo == 1 % if we have to reprocess
        if varOpen.ave.qtype<0, % case 1: use behavior trigger, from roi = abs(varOpen.ave.qtype) index
        bei = seed_info.be(-varOpen.ave.qtype,:);
        
        if varOpen.ave.temp_filter~=0 % if temporal filter present
            beif = smooth(bei,varOpen.ave.temp_filter*length(seed_info.bxx)/max(seed_info.bxx)); % applying filter to data and assigning it to beif 
        else, 
            beif = bei; % temporal filter not present, assign data to beif without smoothing   
        end

        m = mean(beif); % mean of behavior
        s = std(beif); % standard deviation of behaviour
        y = m; % default starting value for clicking = mean
        x = 1; 
        
        %subplot(111) % creates 
        subtightplot(1,1,1,[0.1,0.1])
        while ((x>0)&(x<max(seed_info.bxx))) % while clicking inside the graph window

            plot(seed_info.bxx,beif,'k') % behavior filtered
            hold on
            plot(seed_info.bxx,m*ones(size(seed_info.bxx)),'b'); % line of mean
            plot(seed_info.bxx,(m+s)*ones(size(seed_info.bxx)),'c'); % line of mean + 1SD
            plot(seed_info.bxx,(m+2*s)*ones(size(seed_info.bxx)),'g'); % line of mean + 2SD
            plot(seed_info.bxx,y*ones(size(seed_info.bxx)),'r'); % current threshold
            text(0,double(m),'mean') % differentiate what each different coloured line represents on the graph
            text(0,double(m+s),'mean+SD')
            text(0,double(m+2*s),'mean+2SD')
            text(0,double(y),'threshold')

            ben = bwlabel(beif>y); % period above threshold
            clear benpos benposdiff 
            for ii = 1:max(ben) % for each of selected threshold period
                beni = ben==ii; 
                benpos(ii) = min(find(beni==1)); % find the time (frame) location of the first frame
            end % and store the frame tag in benpos

            for ii = 2:length(benpos) % for each of the frame tag, calculate the interval
                benposdiff(ii) = (benpos(ii)-benpos(ii-1))*max(seed_info.bxx)/length(seed_info.bxx);
            end
            benposdiff(1) = Inf; % the first interval in not known
            benpos = benpos(find(benposdiff>varOpen.ave.refr_period)); % only keep the tag with interval larger than refractory period

            for ii = 1:length(benpos) % for each tag, place it on the graph
                plot(benpos(ii)*max(seed_info.bxx)/length(seed_info.bxx),y,'ro')
            end

            hold off
            disp('CLICK ON THRESHOLD AND THEN OUTSIDE')

            [x,y]=ginput(1);
        end
        tps = benpos*max(seed_info.bxx)/length(seed_info.bxx); % convert in real time (s)

    else, % case 2: use text file

        switch varOpen.ave.qtype % if positive, it will locate go trials through using textfiles 
            case 2, qstring = 'GO=2'; % look for 'GO=2' in text
            case 3, qstring = 'GO=-4'; 
            case 4, qstring = 'GO=1';
            case 5, qstring = 'GO=-1';
            case 6, qstring = 'GO=-3';
            case 7, qstring = 'GO=-2'; 
            otherwise, qstring = 'noanychancethisstringcanexist';
        end
        tps = []; % tps is the vector of time to average 
        sum_tps = 0 ; % determines the trial (to be processed) number corresponding to times of occurrence
        for i = 1:length(run_list)
            varOpen.ave.tps(run_list(i)).tps_n = []; % determines the times of occurrence of the qtype chosen
            varOpen.ave.tps(run_list(i)).pertrial = []; % trial number 
            tps_pertrial = [];
            tps_expe=[]; 
            if varOpen.ave.qtype == 0 % focus on data whenever the mouse licked 
                varOpen.ave.tps(run_list(i)).lick = varOpen.varProc.events(run_list(i)).event_lick;
            elseif varOpen.ave.qtype == 1 % focus on data where the mouse was rewarded 
                varOpen.ave.tps(run_list(i)).reward = varOpen.varProc.events(run_list(i)).event_reward;
            elseif varOpen.ave.qtype > 1 % refers to the different GO trials where the mouse responded to a stimulus in a specific way
                    for j = 1:length(varOpen.varProc.events(run_list(i)).event_task_info)
                        if length(findstr(char(varOpen.varProc.events(run_list(i)).event_task_info(j)),qstring)) ~= 0 
                            tps_pertrial = [tps_pertrial varOpen.varProc.events(run_list(i)).event_task(j)] ; % times of GO trial
                            sum_tps = sum_tps +1;
                            tps_expe= [tps_expe sum_tps] ;
    %                         varOpen.ave.tps_expe = [tps_expe varOpen.varProc.run_list(varOpen.varProc.run_list(i))] ;
                        end 
                    end 
                    if i == length(run_list)
                        for j = 1:length(varOpen.varOutput.event_list_task_info)
                            if length(findstr(char(varOpen.varOutput.event_list_task_info(j)),qstring))~= 0 % find all events with qtype specified 
                                tps = [tps varOpen.varOutput.event_list_task(j)]; % time of occurrence for specified qtype 
                            end
                        end
                    end
                    if length(tps_pertrial) == 0 
                        continue 
                    else 
                        varOpen.ave.tps(run_list(i)).pertrial = tps_pertrial; % times of GO trial for every expe 
                        varOpen.ave.tps(run_list(i)).tps_n= tps_expe; % The trial n each cue corresponds to 
                    end
            end
            end
        end
        
        varOpen.Output.tps = tps;
    % exclude trials that exceed the user defined threshold of the specified ROI

    if varOpen.ave.mvt_roi ~= 0 % if mvt discarder active
        bei = seed_info.be(varOpen.ave.mvt_roi,:); % behavioural data of specified ROI
        if varOpen.ave.temp_filter~=0 % if filter smoothing present
            beif = smooth(bei,varOpen.ave.temp_filter*length(seed_info.bxx)/max(seed_info.bxx));
        else, beif = bei; end
        m = mean(beif); % mean of behavior
        s = std(beif); % standard deviation of it   
        seuil = m + varOpen.ave.mvt_kSD*s; % seuil = threshold determined by adding kSD times of std dev to mean
        mvt_tag = beif>seuil; % place 1 where behaviour exceeds threshold 
        
        for i = 1:length(run_list) 
            bei_expe = varOpen.Seed_expe(run_list(i)).be_expe(varOpen.ave.mvt_roi,:); % take behavioural data of the specified ROI 
            if varOpen.ave.temp_filter ~= 0 
                beif_expe = smooth(bei_expe,varOpen.ave.temp_filter*length(seed_info.bxx)/max(seed_info.bxx)); % filter the behavioural data
            else 
                beif_expe = bei ; 
            end
            varOpen.ave.expe(run_list(i)).mvt_tag_expe = beif_expe >seuil ; % have matrix stating whether or not the behavioural data exceeds the threshold 
        end 
                
        figure(2) % Indicates the times in which the behavioural data exceeds the threshold as well as the time each cue was given 
        plot(seed_info.bxx,beif,'k'), hold on
        plot(seed_info.bxx,.9*mvt_tag*max(beif),'r'),
        for i = 1:length(tps)
            plot(tps(i),max(beif),'vb');
            title('Calcium Spiking In Relation To Its Movement')
            xlabel('Time(s)')
        end
        hold off
        
    else
        mvt_tag = zeros(sum(varOpen.varOutput.jonction),1); % else, vector of zeros
        for i = 1:length(run_list) % if filtering isn't needed, take into account all trials (no threshold) 
            varOpen.ave.expe(run_list(i)).mvt_tag_expe = zeros(size(varOpen.Seed_expe(run_list(i)).be_expe,2),1) ; 
        end 
    end
   
    I2 = []; % initiate the matrix of averaging [X,Y,Frame,Channel,Trial]
    I2_be = [];% initiate the matrix of averaging of behavior mov [Frame,Trial,roi]
    lick_count = zeros(sum_tps,(varOpen.ave.pre/0.2)*2-1); % takes into account licks happening every 0.2s for each expe and gets avged 
    not_added = []; % trials that should not be added 
    
        for i = 1:length(run_list)
            preSF_expe= round(expe(run_list(i)).realSF/varOpen.varProc.tb*varOpen.ave.pre); % determine n frames before and after cue 
            preSF_beh = round(expe(run_list(i)).realSF_cam/varOpen.varProc.tb*varOpen.ave.pre); % both calcium and behavioural data have different sampling frequencies, these values could be different
            if length(varOpen.ave.tps(run_list(i)).pertrial) == 0 % if there are no cues of specified qtype in the run_list, proceed to next and skip 
                continue 
            else 
            for j = 1:length(varOpen.ave.tps(run_list(i)).pertrial)
                lick_tps = varOpen.varProc.events(run_list(i)).event_lick; % times of every lick per trial
                pos_expe = round(varOpen.ave.tps(run_list(i)).pertrial(j)*(expe(run_list(i)).realSF/varOpen.varProc.tb));% position of fluo frame at time specified in per trial
                pos_beh = round(varOpen.ave.tps(run_list(i)).pertrial(j)*(expe(run_list(i)).realSF_cam/varOpen.varProc.tb)); % position of beh frame at specified time 

            try
                I = expe(run_list(i)).If(:,:,pos_expe-preSF_expe:pos_expe+preSF_expe,:); % cut calcium data of trial from expe
                mvt_tag_cut = varOpen.ave.expe(run_list(i)).mvt_tag_expe(pos_beh-preSF_beh:pos_beh+preSF_beh); % indicates frames that exceed the threshold with "1"
                mvt_tag_cut = mvt_tag_cut(1:round(varOpen.ave.mvt_pre*length(mvt_tag_cut))); % cut out duration of frames where we want to see whether or not the frames in that time course exceed threshold. 

                try
                    I_be = varOpen.Seed_expe(run_list(i)).be_expe(:,pos_beh-preSF_beh:pos_beh+preSF_beh); % cut out behaviour data corresponding to the trial being looked at 
                end


                if sum(mvt_tag_cut)==0 % if threshold is exceeded, exclude the trial in averaging
                    for k = 1:size(I,4) % for each channel
                        BL = mean(I(:,:,1:preSF_expe,k),3); % baseline calculation
                        for f = 1:size(I,3) % for each frame
                            I(:,:,f,k) = I(:,:,f,k)-BL; % subtract the baseline
                            % NB: I is already DF/F so no need to devide again
                            %     substraction was mostly to remove eventually
                            %     fluctuation of baseline
                        end
                    end
                    
                    if length(I2) ~= 0 % If there already are concatenated calcium elements 
                        if size(I,3) ~= size(I2,3) % If frames (time) of calcium data aren't the same, reshape it so number of frames and times are consistent 
                            for f = 1:size(I,4)
                                I_reshape = I(:,:,:,f); 
                                if f == 1
                                I_new = [];
                                else
                                end 
                                I_reshape = reshape(I_reshape,[size(I_reshape,1)*size(I_reshape,2), size(I_reshape,3)]);
                                I_reshape = imresize(I_reshape,[size(I_reshape,1), size(I2,3)]); 
                                I_reshape = reshape(I_reshape,[size(I2,1),size(I2,2),size(I2,3)]); 
                                I_new = cat(4,I_new,I_reshape); 
                            end 
                            I = I_new;
                        else 
                        end 
                    else 
                    end 
                    I2 = cat(5,I2,I); % concatenate the trials by the 5th dimension, indicating n trials that meet this criteria 
                    for c = 0:((varOpen.ave.pre/0.2)*2-1) % want 0.2s increments for n of licks for both prestim and poststim periods so it can be averaged (show whether or not it aligns with calcium and movement data) 
                        n_licks = find(((varOpen.ave.tps(run_list(i)).pertrial(j)-varOpen.ave.pre+0.2*c)< lick_tps) & (lick_tps < (varOpen.ave.tps(run_list(i)).pertrial(j)-varOpen.ave.pre+0.2*(c+1))));
                        lick_count(varOpen.ave.tps(run_list(i)).tps_n(j),c+1) = length(n_licks)  ; % take into account licks happening in every expe every 0.2s, put into matrix  
                    end

                    try
                        I2_be = cat(3,I2_be,I_be); % concatenate the behavioural data 
                        
                    end

                    disp([num2str(size(I2,5)) ' trials added']) 
                else 
                    disp('trial not added (because between 2 trials)')
                    not_added = [not_added varOpen.ave.tps(run_list(i)).tps_n(j)] ;
                end
            catch
                disp('trial cannot be added')
                not_added= [not_added varOpen.ave.tps(run_list(i)).tps_n(j)] ;

            end
            end 
            end
        end
    lick_count(not_added,:) = []; % exclude licks of the trials that could not be added 
varOpen.ave.lick_count= lick_count ; % save all licks as a field 
    
    disp(['number of trials: ' num2str(size(I2,5))]) 
    disp(['qtype=' num2str(varOpen.ave.qtype)])
    disp('   0: lick')
    disp('   1: reward')
    disp('   2: GO=2')
    disp('   3: GO=-4')
    disp('   4: GO=1')
    disp('   5: GO=-1')
    disp('   6: GO=-3')
    disp('   7: GO=-2')
    
    disp(['saving ' Ave_temp_file])
    if length(I2_be~=0) % save all data in ave_temp so it can be accessed again without having to reprocess everything (with the same criteria) 
        save(Ave_temp_file,'I2','I2_be','lick_count','qstring'); % save average in temp file, with behavior roi if here
        save(varOpen.ave.outfile,'I2','I2_be','lick_count','qstring')
    else
        save(Ave_temp_file,'I2','qstring'); % save average in temp file
        save(varOpen.ave.outfile,'I2','qstring')
    end
else 
end

    if isempty(varOpen.ave.ntrials)== 1   % If range of trials aren't defined, use the whole range. 
    varOpen.ave.ntrials= [1:size(I2,5)]; %I2 specifies x,y,n-frames, channel, trials 
    else 
    end 

    try % only include the trials specified in the master script (if not the full sequence) 
        I3 = I2(:,:,:,1,varOpen.ave.ntrials) - I2(:,:,:,2,varOpen.ave.ntrials); % eventually substract the reflectance
    catch
        I3 = I2;
    end
    I3 = squeeze(I3); % squeeze the 4th dimension who doesnt exist anymore
    I4 = mean(I3,4); % average trial
    
    I4 = OIA_multiDC(I4, varOpen.varOutput.roi); % add roi to the average between trial
    
    %I4(1:10,1:10,preSF+1) = max(I4(:)).*ones(10); % add a white square at the go-trial
    
    I4(1,1,:)=max(I4(:)).*ones(1,1,size(I4,3)); % for normalization, put min and max in left-up corner
    I4(1,2,:)=min(I4(:)).*ones(1,1,size(I4,3));
    

clear I5 % creation of the montage of the response

for i = 1:10
  I5(:,:,i) = max(I4(:,:,round((i-1)*size(I4,3)/10)+1:round(i*size(I4,3)/10)),[],3);
  %I5(:,:,i) = mean(I4(:,:,round((i-1)*size(I4,3)/10)+1:round(i*size(I4,3)/10)),3);
end

MON = [];
if length(varOpen.ave.range)==0
    montage_range = [min(I5(:)) max(I5(:))];
else
    montage_range = [varOpen.ave.range];
end
for i = 1:size(I5,3) % for each frame of the I5
    
    MONi = I5(:,:,i);
    MONi(1,1) = montage_range(1);
    MONi(1,2) = montage_range(2);
    
    
    MONi = OIA_color(MONi,jet); % convert in color
    MONi = OIA_multiDC(MONi,uint8(varOpen.varOutput.roi)); % apply the mask
    MON = cat(2,MON,uint8(MONi)); % concatenate the results in x-axis
end

figure(1)
fig1 = figure(1);
fig1.Renderer = 'Painters'; % vector graphics 
fig1.Position = [680 505 784 473]; % specifying dimension the figure should be displayed as 
% subplot(311), imshow(MON) % and display the montage
subtightplot(3,2,[1 2],[0.01,0.01]), imshow(MON) % subtightplot not a built in matlab function; subtightplot.m can be found in OIA 
txx_montage = [-varOpen.ave.pre:varOpen.ave.pre/5:varOpen.ave.pre];

for i = 1:11
    text(size(I4,1)*(i-1),-size(I4,1)/10,num2str(txx_montage(i)), 'FontSize', 13)
end
colormap jet
title({'Montage of Brain Images Corresponding to the Time (in Seconds)'; ' ' }, 'FontSize', 15)

map = max(I4(:,:,26:end),[],3)-min(I4(:,:,1:25),[],3); % taking the range of the time dimension (frame) 
disp(size(I4))
%map = max(I4,[],3)-min(I4,[],3); % taking the range of the time dimension (frame)
if length(varOpen.ave.seed_pixel) == 2 % for specified pixel 
    x= varOpen.ave.seed_pixel(1); % x and y coordinates of seed pixel 
    y= varOpen.ave.seed_pixel(2); 
    subtightplot(3,2,[3 5],[0.01,0.01])
    if length(varOpen.ave.range) == 0
        imshow(map,[min(map(:)) max(map(:))]); % show the max-min of DF/F% 
    else
        
        try
            imshow(map,[varOpen.ave.range(1) varOpen.ave.range(2)]); % show specified range of DF/F %
        catch
             imshow(map,[min(map(:)) max(map(:))]) % show DF/F% over the max and min 
        end
    end
    colormap jet, 
    colorbar
    hold on
    plot(x,y,'ow'); hold off % display the pixel 
    title(['Trials: ' num2str(min(varOpen.ave.ntrials)) '-' num2str(max(varOpen.ave.ntrials)) ' For Pixel ' num2str(x) 'x' num2str(y)], 'FontSize', 15) 
        

    
    s_fluo = squeeze(I2(y,x,:,1,varOpen.ave.ntrials)); % fluo data
    s_fluo_m = mean(s_fluo,2);
    s_fluo_s = std(s_fluo,[],2)./sqrt(size(s_fluo,2));
    try,
        s_refl = squeeze(I2(y,x,:,2,varOpen.ave.ntrials)); % refl data
        s_refl_m = mean(s_refl,2);
        s_refl_s = std(s_refl,[],2)./sqrt(size(s_refl,2));
    end
    
    txx = -1+2*[1:length(s_fluo_m)]'/length(s_fluo_m); % time points for the calcium data
    txx = txx*varOpen.ave.pre; % taking into account the time we want to display for pre and post stimulus 
    
    subtightplot(3,2,6,[0.05,0.05]) % display mean and m +- standard deviation 
    plot(txx,s_fluo_m,'g','LineWidth',2); hold on % fluo
    plot(txx,s_fluo_m-s_fluo_s,'g'); 
    plot(txx,s_fluo_m+s_fluo_s,'g'); 
    try
        plot(txx,s_refl_m,'b','LineWidth',2) % refl
        plot(txx,s_refl_m-s_refl_s,'b')
        plot(txx,s_refl_m+s_refl_s,'b')
    end

        lick_count = mean(lick_count(varOpen.ave.ntrials,:)); % average all licks for each 0.2s interval to graph 
        lick_count = lick_count./0.2 ;
        lick_xxplot= [] ; 
        lick_xplot = -varOpen.ave.pre:0.2: varOpen.ave.pre; % take the average of licks between all trials and plot it
            for t = 1:length(lick_xplot)-1 
                lick_xxplot(t) = (lick_xplot(t) + lick_xplot(t+1))/2;
            end
            
try
            yyaxis right 
            ylabel('licks/sec')  %left axis represents DF/F% data and right represents the axis for licks 
            plot(lick_xxplot, lick_count) 
    end 

    if varOpen.ave.subtraction_line == 1 % states whether or not to display fluo-refl data
 
        sub_mean= s_fluo_m-s_refl_m; % corrected mean data 
        sub_max= (s_fluo_m-s_fluo_s)-(s_refl_m-s_refl_s); % mean- stddev of corrected signals  
        sub_low= (s_fluo_m+s_fluo_s)-(s_refl_m+s_refl_s); % mean + stddev of corrected signals 
        
    try 
        yyaxis left 
        plot(txx,sub_mean,'k','LineWidth',2)
        plot(txx,sub_max, 'k')
        plot(txx,sub_low, 'k')
    end 
    else 
    end
  
    for i = 1:length(txx_montage) % every line indicates the corresponding time corresponding to the montage shown at the very top
        plot([txx_montage(i) txx_montage(i)],[min(s_fluo_m) max(s_fluo_m)],':r','LineWidth',.5)
    end
       

    try
        subtightplot(3,2,4,[0.05,0.05])
        clear s_be_m_liste % plot all ROI movements for the specified range of trials 
        for i = 1:size(I2_be,1)
            s_be = squeeze(I2_be(i,:,varOpen.ave.ntrials));
            s_be_m_liste(i) = max(mean(s_be,2));
        end
        s_be_m_liste = max(s_be_m_liste);
        for i = 1:size(I2_be,1) % concerning movement ROIs 
            s_be = squeeze(I2_be(i,:,varOpen.ave.ntrials)); 
            s_be_m = mean(s_be,2); % mean 
            s_be_s = std(s_be,[],2)/sqrt(size(s_be,2)); % standard deviation 
            % s_be_max = max(s_be_m);
%             s_be_m = (max(s_fluo_m).*s_be_m/s_be_max);
%             s_be_s = max(s_fluo_s).*s_be_s/s_be_max;
            % This now normalize by the maximum movement value quantified
            % for all roi: s_be_m_liste
            s_be_m = (max(s_fluo_m).*s_be_m/s_be_m_liste);
            s_be_s = max(s_fluo_s).*s_be_s/s_be_m_liste;
            
            bxx= varOpen.ave.pre*(-1+2*[1:length(s_be_m)]'/length(s_be_m));
            
            s_be_m_i = interp1(bxx,s_be_m,txx); % interpolate so bxx has same time points as txx 
            s_be_m_i(isnan(s_be_m_i))= 0; 
            s_be_m_l = interp1(bxx,s_be_m-s_be_s,txx); 
            s_be_m_l(isnan(s_be_m_l))= 0; 
            s_be_m_h = interp1(bxx,s_be_m+s_be_s,txx); 
            s_be_m_h(isnan(s_be_m_h))= 0; 
            
            plot(txx,s_be_m_i,'Color',seed_info.couleur(i,:),'LineWidth',2)
            hold on
            plot(txx,s_be_m_l,'Color',seed_info.couleur(i,:)) % plot the stddev as well as the mean for each ROI
            plot(txx,s_be_m_h,'Color',seed_info.couleur(i,:)) 
            
            text(double(txx(end)),double(s_be_m(end)),['roi' num2str(i)],'Color',seed_info.couleur(i,:))

        end
    try, title(qstring); end    
    end
    
    
    hold off
    subtightplot(3,2,6,[0.05,0.05])
    yyaxis left 
    xlabel('time(s)')
    ylabel('DF/F')
    
    if length(varOpen.ave.range) == 0
        range = [nanmax(s_fluo_m) - nanmin(s_fluo_m)];
        varOpen.ave.range(1) = nanmin(s_fluo_m) - .5*range;
        varOpen.ave.range(2) = nanmax(s_fluo_m) + .5*range;
    end
     
    try,
        ylim([varOpen.ave.range(1) varOpen.ave.range(2)])
    end
    
else 
map=max(I4, [],3)- min(I4,[],3); 
x=1; y=1; button = 1; % for the interactive map; defined as first pixel initially, then myginput is executed and x y will change 

while button ~= 3 
    subtightplot(3,2,[3 5],[0.01,0.01])
    if length(varOpen.ave.range) == 0
        imshow(map,[min(map(:)) max(map(:))]); 
    else
        
        try
            imshow(map,[varOpen.ave.range(1) varOpen.ave.range(2)]); 
        catch
            imshow(map,[min(map(:)) max(map(:))]); 
        end
    end
    colormap jet, 
    colorbar
    hold on
    plot(x,y,'ow'); hold off
    title(['Trials: ' num2str(min(varOpen.ave.ntrials)) '-' num2str(max(varOpen.ave.ntrials))], 'FontSize', 15)
    
    s_fluo = squeeze(I2(y,x,:,1,varOpen.ave.ntrials)); % get the fluo data
    s_fluo_m = mean(s_fluo,2);
    s_fluo_s = std(s_fluo,[],2)./sqrt(size(s_fluo,2));
    try,
        s_refl = squeeze(I2(y,x,:,2,varOpen.ave.ntrials));
        s_refl_m = mean(s_refl,2);
        s_refl_s = std(s_refl,[],2)./sqrt(size(s_refl,2));
    end
    
    txx = -1+2*[1:length(s_fluo_m)]'/length(s_fluo_m);
    txx = txx*varOpen.ave.pre;
    subtightplot(3,2,6,[0.05,0.05])
    plot(txx,s_fluo_m,'-g','LineWidth',2); hold on
    plot(txx,s_fluo_m-s_fluo_s,'-g'); 
    plot(txx,s_fluo_m+s_fluo_s,'-g'); 
    try
        plot(txx,s_refl_m,'-b','LineWidth',2)
        plot(txx,s_refl_m-s_refl_s,'-b')
        plot(txx,s_refl_m+s_refl_s,'-b')
    end
    
            lick_count_m = mean(lick_count(varOpen.ave.ntrials,:)); 
        lick_count_m = lick_count_m./0.2 ;
        lick_xxplot= [] ; 
        lick_xplot = -varOpen.ave.pre:0.2: varOpen.ave.pre; % this killed it. 
            for t = 1:length(lick_xplot)-1 
                lick_xxplot(t) = (lick_xplot(t) + lick_xplot(t+1))/2;
            end
    
    if varOpen.ave.subtraction_line == 1
        
    try 
        plot(txx,s_fluo_m-s_refl_m,'-k','LineWidth',2)
        plot(txx,(s_fluo_m-s_fluo_s)-(s_refl_m-s_refl_s), '-k')
        plot(txx,(s_fluo_m+s_fluo_s)-(s_refl_m+s_refl_s), '-k')
    end

    else 
    end
    
    for i = 1:length(txx_montage)
        plot([txx_montage(i) txx_montage(i)],[min(s_fluo_m) max(s_fluo_m)],':r','LineWidth',.5)
    end
hold off 
    try
            yyaxis right 
            ylabel('licks/sec') 
            plot(lick_xxplot, lick_count_m, '-') 
            yyaxis left 
    end
    
    try
        subtightplot(3,2,4,[0.05,0.05])
        clear s_be_m_liste
        for i = 1:size(I2_be,1)
            s_be = squeeze(I2_be(i,:,varOpen.ave.ntrials));
            s_be_m_liste(i) = max(mean(s_be,2));
        end
        s_be_m_liste = max(s_be_m_liste);
        for i = 1:size(I2_be,1)
            s_be = squeeze(I2_be(i,:,varOpen.ave.ntrials)); % behavioural data for specified trials 
            s_be_m = mean(s_be,2);
            s_be_s = std(s_be,[],2)/sqrt(size(s_be,2));
            % s_be_max = max(s_be_m);
%             s_be_m = (max(s_fluo_m).*s_be_m/s_be_max);
%             s_be_s = max(s_fluo_s).*s_be_s/s_be_max;
            % This now normalize by the maximum movement value quantified
            % for all roi: s_be_m_liste
            s_be_m = (max(s_fluo_m).*s_be_m/s_be_m_liste);
            s_be_s = max(s_fluo_s).*s_be_s/s_be_m_liste;
            
            bxx= (-1+2*[1:length(s_be_m)]'/length(s_be_m))*varOpen.ave.pre;

            s_be_m_i = interp1(bxx,s_be_m,txx); 
            s_be_m_i(isnan(s_be_m_i))= 0; 
            s_be_m_l = interp1(bxx,s_be_m-s_be_s,txx); 
            s_be_m_l(isnan(s_be_m_l))= 0; 
            s_be_m_h = interp1(bxx,s_be_m+s_be_s,txx); 
            s_be_m_h(isnan(s_be_m_h))= 0; 
            
            plot(txx,s_be_m_i,'Color',seed_info.couleur(i,:),'LineWidth',2)
            hold on
            plot(txx,s_be_m_1,'Color',seed_info.couleur(i,:)) % plot the stddev as well as the mean for each ROI
            plot(txx,s_be_m_h,'Color',seed_info.couleur(i,:)) 
            
            text(double(txx(end)),double(s_be_m(end)),['roi' num2str(i)],'Color',seed_info.couleur(i,:))
            
        end
       try title(qstring); end 
    end
    
    hold off
    subtightplot(3,2,6,[0.05,0.05])
    yyaxis left
    xlabel('time(s)')
    ylabel('DF/F')
%     hold off 
    
    
    if length(varOpen.ave.range) == 0
        range = [nanmax(s_fluo_m) - nanmin(s_fluo_m)];
        varOpen.ave.range(1) = nanmin(s_fluo_m) - .5*range;
        varOpen.ave.range(2) = nanmax(s_fluo_m) + .5*range;
    end

    try,
        ylim([varOpen.ave.range(1) varOpen.ave.range(2)])
    end

    [x,y,button]=myginput
    while x>size(I2,1), x =x-size(I2,1); end 

end
varOpen.ave.tps_expe_trials = tps_expe_trials;
end



