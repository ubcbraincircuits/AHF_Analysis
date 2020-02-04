% Main analysis pipeline for figures 9 & 10
% This code works on the output from the FARM pipeline
%
%
clear, clc

%% path information
% directories containing data
path_info.go_hit          = 'B:\AHF\Code+2\';
path_info.go_late         = 'B:\AHF\Code-2\';
path_info.go_early        = 'B:\AHF\Code-4\';

% directories where results will be saved
path_info.save_path           = 'B:\AHF\Outputs\';

%% load files

go_hit_files    =  helper.getAllFiles( path_info.go_hit );
go_late_files   =  helper.getAllFiles( path_info.go_late );
go_early_files  =  helper.getAllFiles( path_info.go_early );

assert(length(go_hit_files) == length(go_late_files) && ...
    length(go_hit_files) == length(go_early_files), ...
    'Number of animals in each condition should be equal.');

num_animals = length(go_hit_files);

% go hit
for i = 1:num_animals
    go_hit(i).data = load([path_info.go_hit, go_hit_files{i}],'I','be_all');
    go_hit(i).ID = go_hit_files{i}(2:end-4);   
    go_hit(i).roi = get_roi_from_data(go_hit(i).data.I(:,:,1,1));
    
    go_late(i).data = load([path_info.go_late, go_late_files{i}],'I','be_all');
    go_late(i).ID = go_late_files{i}(2:end-4);
    go_late(i).roi = get_roi_from_data(go_late(i).data.I(:,:,1,1));
    
    go_early(i).data = load([path_info.go_early, go_early_files{i}],'I','be_all');
    go_early(i).ID = go_early_files{i}(2:end-4);
    go_early(i).roi = get_roi_from_data(go_early(i).data.I(:,:,1,1));
end
disp('Finished loading data')

%% register ROIs

% --------------------- choose a reference image -------------------------
plotfigs = [];
imsz = size(go_hit(1).roi,1);

for i = 1:num_animals
    temp_plot = cat(1, go_hit(i).roi, go_early(i).roi, go_late(i).roi);        
    plotfigs = cat(2,plotfigs,temp_plot);   
end

figure, imagesc(plotfigs), colormap gray
title('CLICK ON REFERENCE IMAGE');
[x,y] = ginput(1); close gcf;
x = round(x); y = round(y);
refnum = ceil(x/imsz);
refgroup = ceil(y/imsz);

switch refgroup
    case 1
        fixed = go_hit(refnum).roi;
    case 2
        fixed = go_early(refnum).roi;
    case 3
        fixed = go_late(refnum).roi;
    otherwise
end

% ------------------------- register images ----------------------------
for i = 1:num_animals
    
    go_hit(i).reg.I = reg_image_by_roi(go_hit(i).roi, fixed, ...
        go_hit(i).data.I);
    
    go_early(i).reg.I = reg_image_by_roi(go_early(i).roi, fixed, ...
        go_early(i).data.I);
    
    go_late(i).reg.I = reg_image_by_roi(go_late(i).roi, fixed, ...
        go_late(i).data.I);  
end
disp('Image registration complete')

%% choose ROIs

ROIS = {'V1','HL','FL','M1','RS','M2','BC','ALM'};
scaling = 8.6/64; % field of view (mm/px)
[CL,CR,bregma] = get_coords(fixed, scaling, ROIS);

num_rois = length(ROIS);
trial_length = size(go_hit(1).reg.I, 3);
radius=2;

trace_array = cell(num_animals, num_rois);

for j = 1:num_animals
    for i = 1:num_rois
        subplot(num_animals, num_rois, i+(j-1)*num_rois), hold on
        
        % average data in square region around ROI center
        regionyL = CL(i,2)-radius:CL(i,2)+radius;
        regionyR = CR(i,2)-radius:CR(i,2)+radius;
        regionxL = CL(i,1)-radius:CL(i,1)+radius;
        regionxR = CR(i,1)-radius:CR(i,1)+radius;
        
        data_in_roi = go_hit(j).reg.I([regionyL,regionyR],[regionxL,regionxR],:,:);
        all_trials = squeeze( mean(mean(data_in_roi,2),1) ); 
        
        % create cell array containing all single trial data
        trace_array{j,i} = all_trials;
    end
end

%% brain-behavior correlations
rng(7913483)
fs = trial_length/5; % frame rate (downsampled from Matthieu's code)

% group animals by cage
nameorder = {'0629','0797','0965','0252','8423','8474','0205','0219','0270','0240','1234','1245','1379'};
trial_IDs = cell(1,num_animals);
for i = 1:num_animals
    trial_IDs{i} = go_hit(i).ID;
end

% get minimum number of trials of any single mouse
minsize = min(min(cellfun('size', trace_array, 2)));

% preallocate
ALM_z = nan(num_animals, trial_length);    tongue_z    = nan(size(ALM_z));       
HL_z  = nan(num_animals, trial_length);    hindlimb_z  = nan(size(ALM_z));     
BC_z  = nan(num_animals, trial_length);    barrel_z    = nan(size(ALM_z));       

ALM_r = nan(1, num_animals);
HL_r  = nan(1, num_animals);
BC_r  = nan(1, num_animals);

random_index = nan(num_animals, minsize);

for i = 1:num_animals    
    a = contains(trial_IDs, nameorder{i});
    random_index(i,:) = randperm(size(trace_array{a,6},2),minsize);

    % behavior traces
    tongue_z(i,:)   = squeeze( mean( zscore( ...
        go_hit(a).data.be_all(1,:,random_index(i,:)), [], 2), 3) );
    barrel_z(i,:)   = squeeze( mean( zscore( ...
        go_hit(a).data.be_all(2,:,random_index(i,:)), [], 2), 3) );
    hindlimb_z(i,:) = squeeze( mean( zscore( ...
        go_hit(a).data.be_all(3,:,random_index(i,:)), [], 2), 3) );
    
    % fluorescence traces
    ALM_z(i,:) = mean(zscore(trace_array{a,6}(:,random_index(i,:))),2);
    HL_z(i,:)  = mean(zscore(trace_array{a,2}(:,random_index(i,:))),2);
    BC_z(i,:)  = mean(zscore(trace_array{a,8}(:,random_index(i,:))),2);

    % matched brain-behavior correlations
    rval = corrcoef(tongue_z(i,:), ALM_z(i,:));   ALM_r(i) = rval(2,1);    
    rval = corrcoef(hindlimb_z(i,:), HL_z(i,:));  HL_r(i) = rval(2,1);    
    rval = corrcoef(barrel_z(i,:), BC_z(i,:));    BC_r(i) = rval(2,1);
end

% scrambled correlations
ALM_r_scrambled = [];
HL_r_scrambled = [];
BC_r_scrambled = [];

for i = 1:size(ALM_z,1)
    for j = 1:size(ALM_z,1)
        if i~= j && i < j
            rval = corrcoef(ALM_z(i,:), tongue_z(j,:));
            ALM_r_scrambled = cat(1, ALM_r_scrambled, rval(2,1));

            rval = corrcoef(HL_z(i,:), hindlimb_z(j,:));
            HL_r_scrambled = cat(1, HL_r_scrambled, rval(2,1));

            rval = corrcoef(BC_z(i,:), barrel_z(j,:));
            BC_r_scrambled = cat(1, BC_r_scrambled, rval(2,1));
        end
    end
end

% populate array with brain-behavior correlation values
brain_behavior_r = nan(1000,6);
brain_behavior_r(1:length(ALM_r),1) = ALM_r;
brain_behavior_r(1:length(ALM_r_scrambled),2) = ALM_r_scrambled;
brain_behavior_r(1:length(HL_r),3) = HL_r;
brain_behavior_r(1:length(HL_r_scrambled),4) = HL_r_scrambled;
brain_behavior_r(1:length(BC_r),5) = BC_r;
brain_behavior_r(1:length(BC_r_scrambled),6) = BC_r_scrambled;


%% correlation to mean trace
mean_trace_HL = mean(HL_z);
mean_trace_ALM = mean(ALM_z);

single_trial_HL = trace_array(:,2);
single_trial_ALM = trace_array(:,6);

% preallocate
mean_corr_HL  = nan(1,num_animals);    sem_corr_HL  = nan(1,num_animals);
mean_corr_ALM = nan(1, num_animals);   sem_corr_ALM = nan(1, num_animals);

k_ALM = nan(1, num_animals);   k_HL = nan(1, num_animals);
p     = nan(1, num_animals);   h = nan(1, num_animals);

for i = 1:num_animals
    a = contains(trial_IDs,nameorder{i});
    
    num_trials = size(single_trial_HL{a}, 2);
    corrs_HL = nan(1, num_trials);
    corrs_ALM = nan(1, num_trials);
        
    for j = 1:size(single_trial_HL{a},2)
        r_val = corrcoef(single_trial_HL{a}(:,j),mean_trace_HL);
        corrs_HL(j) = r_val(2,1);
        
        r_val = corrcoef(single_trial_ALM{a}(:,j),mean_trace_ALM);
        corrs_ALM(j) = r_val(2,1);
    end
    mean_corr_HL(i) = mean(corrs_HL);
    sem_corr_HL(i) = sem(corrs_HL,2);
    
    
    mean_corr_ALM(i) = mean(corrs_ALM);
    sem_corr_ALM(i) = sem(corrs_ALM,2);
    
    % stats
    k_ALM(i) = kstest(zscore(corrs_ALM));
    k_HL(i) = kstest(zscore(corrs_HL));
    [p(i),h(i)] = signrank(corrs_HL,corrs_ALM);
end

%% correlation figure

% ------------------------------ panel A -------------------------------
figure('Renderer', 'painters', 'Position', [100 300 650 400]),
subplot(2,5,1)
for i = 1:size(HL_z)
    plot(xt(HL_z,fs)-2.5, HL_z(i,:)-i,'k'), hold on
    axis([-3 2.5 -15 0.5]), axis off, title('HL'), 
    line([-2.8 -2.8], [-14.2 -13.2], 'color', [0 0 0]),
    line([-2.8 -1.8], [-14.2 -14.2], 'color', [0 0 0]),
    text(-3.5,-14.4,'1\sigma','Rotation',90,'FontSize',10,'FontWeight','normal');
    text(-3, -14.9, '1 s','FontSize',10,'FontWeight','normal');
end
text(-4.3,2.5,'A','FontSize',14)
subplot(2,5,2)
for i = 1:size(ALM_z)
    plot(xt(ALM_z,fs)-2.5, ALM_z(i,:)-i,'k'), hold on
    axis([-2.5 2.5 -15 0.5]),  axis off, title('ALM')
end

% ----------------------------- panel B ------------------------------
subplot(2,5,6), hold on,
    plot(xt(HL_z,fs)-2.5, mean(HL_z), 'k') 
    shadedErrorBar(xt(HL_z,fs)-2.5, mean(HL_z), sem(HL_z), 'k', 1)
    plot(xt(HL_z,fs)-2.5, mean(hindlimb_z),'b'); 
    shadedErrorBar(xt(HL_z,fs)-2.5, mean(hindlimb_z), sem(hindlimb_z), 'b', 1)
    text(-2,1.05,'Behavior','FontSize',9,'color',[0 0 1])
    text(-2,0.8,'\DeltaF/F_0','FontSize',9,'color',[0 0 0])
    title('HL'), axis([-2.5 2.5 -1 1.2])
    ylabel('Z-score')
    text(-3.7,1.5,'B','FontSize',14)

subplot(2,5,7), hold on,
    plot(xt(HL_z,fs)-2.5,mean(ALM_z),'k') 
    shadedErrorBar(xt(HL_z,fs)-2.5,mean(ALM_z),sem(ALM_z),'k',1)
    plot(xt(HL_z,fs)-2.5,mean(tongue_z),'b'), 
    shadedErrorBar(xt(HL_z,fs)-2.5,mean(tongue_z),sem(tongue_z),'b',1)
    title('ALM'), axis([-2.5 2.5 -1 1.2])
    xlabel('Time (s)')

subplot(2,5,8), hold on,
    plot(xt(HL_z,fs)-2.5,mean(BC_z),'k')
    shadedErrorBar(xt(HL_z,fs)-2.5,mean(BC_z),sem(BC_z),'k',1)
    plot(xt(HL_z,fs)-2.5,mean(barrel_z),'b'), 
    shadedErrorBar(xt(HL_z,fs)-2.5,mean(barrel_z),sem(barrel_z),'b',1)
    title('BC'), axis([-2.5 2.5 -1 1.2])

% ---------------------------- panel C --------------------------------
subplot(2,5,3.2:5.2), hold on

% color background to group mice by cage
    xx = [0 3.5 3.5 0];
    yy = [-0.3 -0.3 1 1];
    patch(xx,yy,'c','FaceAlpha',0.25,'EdgeAlpha',0)

    xx = [3.5 6.5 6.5 3.5];
    yy = [-0.3 -0.3 1 1];
    patch(xx,yy,'m','FaceAlpha',0.25,'EdgeAlpha',0)

    xx = [6.5 9.5 9.5 6.5];
    yy = [-0.3 -0.3 1 1];
    patch(xx,yy,'y','FaceAlpha',0.25,'EdgeAlpha',0)

    xx = [9.5 15 15 9.5];
    yy = [-0.3 -0.3 1 1];
    patch(xx,yy,'k','FaceAlpha',0.15,'EdgeAlpha',0)

h1=errorbar([], mean_corr_HL,  sem_corr_HL,  'bo', 'LineWidth', 1); 
h1.MarkerSize = 3;
h2=errorbar([], mean_corr_ALM, sem_corr_ALM, 'ko', 'LineWidth', 1); 
h2.MarkerSize = 3;

% stats results
h=double(h);
for i = 1:length(h)
    if p(i) < 0.05
        text(i-0.2, double(h(i)*0.9), '*', 'FontSize', 14)
    end
end

line([0 14], [0 0],     'color', [0 0 0], 'LineStyle', ':')
line([0 14], [0.5 0.5], 'color', [0 0 0], 'LineStyle', ':')
xlabel('mouse'), axis([0.5 14 -0.3 1]); 
ylabel('r(single-trial, mean trace)')
legend([h1,h2],{'HL','ALM'},'location','southeast');
text(-1,1.17,'C','FontSize',14);

% ---------------------------- panel D ------------------------------
subplot(2,5,9.5:10), hold on
% arrange matched and mismatched data so violins can be color coded
    matches = brain_behavior_r;
    matches(:,[2,4,6]) = matches(:,[2,4,6])-1000;
    matches = cat(2,matches(:,1:2),1000*ones(size(matches,1),1),matches(:,3:4),...
        1000*ones(size(matches,1),1),matches(:,5:6));
    mismatches = brain_behavior_r;
    mismatches(:,[1,3,5]) = mismatches(:,[1,3,5])-1000;
    mismatches = cat(2,mismatches(:,1:2),1000*ones(size(matches,1),1),mismatches(:,3:4),...
        1000*ones(size(matches,1),1),mismatches(:,5:6));
violinplot(matches,{'ALM-mouth','ALM-mouth(s)','','HL-leg','HL-leg(s)'},'ViolinColor',[0 0 1]);
violinplot(mismatches,{'ALM-mouth','ALM-mouth(s)','','HL-leg','HL-leg(s)'},'ViolinColor',[0 0 0]);
axis([0.25 8.75 -1 1.1])
xticks([1.5,4.5,7.5])
xticklabels({'ALM-tongue','HL-leg','BC-face'})
xtickangle(30)
ylabel('r(\DeltaF/F_0, B)')

% stats
line([1 2], [1.08 1.08], 'color', [0 0 0]);
line([4 5], [1.08 1.08], 'color', [0 0 0]);
line([7 8], [1.08 1.08], 'color', [0 0 0]);
p=signrank(brain_behavior_r(:,1),brain_behavior_r(:,2));
if p >= 0.05
    text(0.5,1.225,'NS')
else
    text(1.15,1.15,'*','FontSize',14)
end
p=signrank(brain_behavior_r(:,3),brain_behavior_r(:,4));
if p >= 0.05
    text(4.15,1.225,'NS')
else
    text(4.15,1.15,'*','FontSize',14)
end
p=signrank(brain_behavior_r(:,5),brain_behavior_r(:,6));
if p >= 0.05
    text(6.5,1.225,'NS')
else
    text(6.5,1.15,'*','FontSize',14)
end

text(-5,1.4,'D','FontSize',14)

%---------------------------- save figure ----------------------------
% saveas(gcf,[path_info.save_path,'correlation_figure2'])

%% ridge regression prediction
rng(121567)

% compile all single trial data from all animals and all ROIs
hit_traces    = cell(1, num_animals);
early_traces  = cell(1, num_animals);
late_traces   = cell(1, num_animals);
for j = 1:num_animals          
    hit_traces{j} = format_data(go_hit(j).reg.I, CL, CR, radius);
    early_traces{j} = format_data(go_early(j).reg.I, CL, CR, radius);
    late_traces{j} = format_data(go_late(j).reg.I, CL, CR, radius);
end
success_trials = zscore(catcell(1,hit_traces),[],3);
early_trials = zscore(catcell(1,early_traces),[],3);
late_trials = zscore(catcell(1,late_traces),[],3);

% compile all dFF maps
single_trial_hit_maps = [];
single_trial_early_maps = [];
single_trial_late_maps = [];

for i = 1:num_animals
    single_trial_hit_maps = cat(4,single_trial_hit_maps,go_hit(i).reg.I);
    single_trial_early_maps = cat(4,single_trial_early_maps,go_early(i).reg.I);
    single_trial_late_maps = cat(4,single_trial_late_maps,go_late(i).reg.I);
end

% averaged single trial maps
average_hit = mean(single_trial_hit_maps,4);
average_early = mean(single_trial_early_maps,4);
average_late = mean(single_trial_late_maps,4);
cat_maps = cat(1,average_hit,average_early,average_late);

% subsample single trial traces and dFF maps
[subsampled_data, index] = subsample_data(success_trials,early_trials,late_trials);

% run decoder for hit vs early
dff_data = catcell(1,subsampled_data{1},subsampled_data{2});
response = [ones(size(dff_data,1)/2,1); -ones(size(dff_data,1)/2,1)];
hitvearly = dffPxlDecoder(dff_data,response);

% run decoder for hit vs late
dff_data = catcell(1,subsampled_data{1},subsampled_data{3});
response = [ones(size(dff_data,1)/2,1); -ones(size(dff_data,1)/2,1)];
hitvlate = dffPxlDecoder(dff_data,response);

%% multiple iterations of full and reduced models using random sampling

rng('default')
num_comparisons = {'hitvearly','hitvlate'};
num_iterations = 1000;
model_types = {'full','-V1','-M2','-ALM','+ALM','+M2','+V1','random'};

for nn = 1:length(num_comparisons)
    % preallocate
    facc     =   nan(trial_length, num_iterations);
    pacc     =   nan(trial_length, num_iterations);
    ALMacc   =   nan(trial_length, num_iterations);
    rALMacc  =   nan(trial_length, num_iterations);
    M2acc    =   nan(trial_length, num_iterations);
    rM2acc   =   nan(trial_length, num_iterations);
    V1acc    =   nan(trial_length, num_iterations);
    rV1acc   =   nan(trial_length, num_iterations);

    for i = 1:num_iterations    
        % subsample single trial traces and dFF maps
        [subsampled_data, index] = subsample_data(success_trials,early_trials,late_trials);
        
        % compare correct vs early and correct vs late
        if nn == 1
            dff_data = catcell(1,subsampled_data{1},subsampled_data{2});
        else
            dff_data = catcell(1,subsampled_data{1},subsampled_data{3});
        end
        response = [ones(size(dff_data,1)/2,1); -ones(size(dff_data,1)/2,1)];

        % random_index for shuffling
        random_index = randperm(size(dff_data,1),size(dff_data,1));

        for mm = 1:length(model_types)
            switch model_types{mm}
                case 'full'
                    idx = [];
                case '-V1'
                    idx = contains(ROIS,'V1');                
                case '-M2'
                    idx = contains(ROIS,'M2');
                case '-ALM'
                    idx = contains(ROIS,'ALM');
                case '+ALM'
                    idx = ~contains(ROIS,'ALM');
                case '+M2'
                    idx = ~contains(ROIS,'M2');
                case '+V1'
                    idx = ~contains(ROIS,'V1');
                case 'random'
                    idx = 1:size(dff_data,2);
            end

            new_dff = dff_data;
            new_dff(:,idx,:) = dff_data(random_index,idx,:);
            dec = dffPxlDecoder(new_dff,response);

            switch model_types{mm}
                case 'full'
                    facc(:,i) = dec.accuracy;
                case '-V1'
                    rV1acc(:,i) = dec.accuracy;
                case '-M2'
                    rM2acc(:,i) = dec.accuracy;
                case '-ALM'
                    rALMacc(:,i) = dec.accuracy;
                case '+ALM'
                    ALMacc(:,i) = dec.accuracy;
                case '+M2'
                    M2acc(:,i) = dec.accuracy;
                case '+V1'
                    V1acc(:,i) = dec.accuracy;
                case 'random'
                    pacc(:,i) = dec.accuracy;
            end

            clear dec
        end
        disp([num2str(i),'/',num2str(num_iterations)])

    end
    
    % save results
%     if nn == 1
%         save([path_info.save_path,'accuracy_results_early.mat'],...
%             'facc','pacc','rALMacc','ALMacc','rM2acc','M2acc','V1acc','rV1acc')
%     else
%         save([path_info.save_path,'accuracy_results_late.mat'],...
%             'facc','pacc','rALMacc','ALMacc','rM2acc','M2acc','V1acc','rV1acc')
%     end
end

%% task decoding figure

figure('Renderer', 'painters', 'Position', [100 300 650 600])
% -------------------------- panel A ------------------------------
map_montage = [];

timepts = 1:5:51;
for i = 1:length(timepts)
    map_montage = cat(2,map_montage,cat_maps(:,:,timepts(i)));
end

subplot(4,4,1:3)
imagesc(map_montage), colormap(gca, 'jet') 
xticks(64*(1:length(timepts))-32), xticklabels(sprintfc('%.1f',timepts/fs - 2.5))
yticks((0:2)*64 + 32), yticklabels({'correct','early','late'})
line([64*find(timepts>=26,1)-64 64*find(timepts>=26,1)-64],[1 size(map_montage,1)],'color',[1 1 1],'LineWidth',2,'LineStyle',':')
line([64*find(timepts>=26,1) 64*find(timepts>=26,1)],[1 size(map_montage,1)],'color',[1 1 1],'LineWidth',2,'LineStyle',':')
line([64*find(timepts>=26,1)-64 64*find(timepts>=26,1)],[size(map_montage,1) size(map_montage,1)],'color',[1 1 1],'LineWidth',2,'LineStyle',':')
line([64*find(timepts>=26,1)-64 64*find(timepts>=26,1)],[1 1],'color',[1 1 1],'LineWidth',2,'LineStyle',':')
set(gca, 'TickLength', [0 0])
xlabel('time (s)')
% c=colorbar; c.Label.String = 'z-scored dF/F';
text(-90,-25,'A','FontSize',14)

% -------------------------- panel B ------------------------------
subplot(4,4,4), imagesc(fixed), hold on, colormap(gca,'gray'), axis off
plot(CL(:,1),CL(:,2),'r*','MarkerSize',2); 
text(CR(:,1)-2,CR(:,2)-0.5,ROIS,'color',[1 0 0],'FontSize',8,'FontWeight','Bold'); 
text(-13, -8,'B','FontSize',14)

% -------------------------- panel C ------------------------------
subplot(4,4,5:7),  hold on
h1=plot(xt(hitvearly.accuracy,fs)-2.5,100*hitvearly.accuracy,'b');
h2=plot(xt(hitvlate.accuracy,fs)-2.5,100*hitvlate.accuracy,'r');
h3=plot(xt(hitvearly.accuracy,fs)-2.5,100*mean(hitvearly.shuffle.accuracy),'color',[0.5 0.5 0.5]);
line([0 0], [45 100], 'color',[0 0 0],'LineStyle',':','LineWidth',1.5)
axis([-2.5 2.5 45 90])
ylabel('accuracy (%)'), xlabel('time (s)'), 
legend([h1,h2,h3],{'correct vs early', 'correct vs late', 'shuffled'},'location','NorthWest','Box','Off')
text(-3.15, 100.5,'C','FontSize',14)

% -------------------------- panel D ------------------------------
subplot(4,4,[8,12,16]), hold on, title('Coefficient Weight')
count=1;
ROI_order = {'ALM','M2','M1','FL','HL','BC','RS','V1'};
for i = 1:length(ROIS)
    line([1 51],[-count*1.2 -count*1.2]+0.2,'color',[0 0 0])
    [~,~,index] = intersect(ROI_order{count}, ROIS);
    h1=plot(hitvearly.weights(:,index)-count*1.2+0.2,'b','LineWidth',1.25);
    h2=plot(hitvlate.weights(:,index)-count*1.2+0.2,'r','LineWidth',1.25);    
    axis([1 53 -10.5 0]), axis off, 

    % trace labels
    text(-3, -((count*1.2)+0.3), ROI_order{count},'Rotation',90,'FontSize',9)
    count = count+1;
end
% scalebar
line([53 53],[-10.1 -9.6], 'color', [0 0 0]), 
line([53 53-fs],[-10.1 -10.1],'color',[0 0 0]),
text(57, -10.1, '0.5', 'Rotation', 90)
text(45, -10.4, '1s')

% cue
line([26 26], [-10 -0.1],'color',[0 0 0],'linestyle',':')

text(-10, 0.63,'D','FontSize',14)

% -------------------------- panel E ------------------------------
subplot(4,4,[9.5:11,13.5:15]), hold on
reduced_model_datasets = {[path_info.save_path,'accuracy_results_early.mat'], ...
    [path_info.save_path,'accuracy_results_late.mat']};

% load reduced model datasets
for i = 1:length(reduced_model_datasets)
    load(reduced_model_datasets{i});
    acc_results = [];
    acc_results(:,1) = max(facc);
    acc_results(:,2) = max(rV1acc);
    acc_results(:,3) = max(rM2acc);
    acc_results(:,4) = max(rALMacc);
    acc_results(:,5) = max(ALMacc);
    acc_results(:,6) = max(M2acc);
    acc_results(:,7) = max(V1acc);
    acc_results(:,8) = max(pacc);
    
    if i == 1
        X = acc_results;
    else
        Y = acc_results;
    end
end

violinplot(X*100,{},'ViolinColor',[0 0 1])
violinplot(Y*100,{},'ViolinColor',[1 0 0])
axis([0 8.75 50 85])
text(0.5,56,'correct vs late','FontSize',8,'Color',[1 0 0])
text(0.5,58,'correct vs early','FontSize',8,'Color',[0 0 1])
labs = {'full','-V1','-M2','-ALM','+ALM','+M2','+V1','random'};
xticklabels(labs), xtickangle(45)
ylabel('Max Accuracy (%)')
text(-4.2,88,'E','FontSize',14)
ax1 = gca;
yyaxis right

ES1 = zeros(1,size(X,2));
ES2 = zeros(1,size(Y,2));
for i = 2:length(labs)
    ES1(i) = calc_effectsize(X(:,1),X(:,i));
    ES2(i) = calc_effectsize(Y(:,1),Y(:,i));
end

plot(ES1,'b')
plot(ES2,'r-')
ax2 = gca;
ax2.YAxis(1).Color = 'k';
ax2.YAxis(2).Color = 'k';
ax2.YAxis(2).Label.String = 'Effect Size';   
ax2.YLim = [0 35];

% stats repeated measures anova
ModelIteration = (1:1000)';
newcell = {'ModelIteration','facc','rV1acc','rM2acc','rALMacc','ALMacc','M2acc','V1acc','random'};

t1 = table(ModelIteration,X(:,1),X(:,2),X(:,3),X(:,4),X(:,5),X(:,6),X(:,7),X(:,8), ...
'VariableNames',newcell);
t2 = table(ModelIteration,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5),Y(:,6),Y(:,7),Y(:,8), ...
'VariableNames',newcell);

ModelType = [1 2 3 4 5 6 7 8]';

rm1 = fitrm(t1,'facc-random ~ ModelIteration','WithinDesign',ModelType);
rm2 = fitrm(t2,'facc-random ~ ModelIteration','WithinDesign',ModelType);

ranovatbl1 = ranova(rm1);
ranovatbl2 = ranova(rm2);

% post hoc bonferroni
[~,~,stats] = anova1(X,[],'off');
c1=multcompare(stats,'CType','bonferroni','Display','off');

[~,~,stats] = anova1(Y,[],'off');
c2=multcompare(stats,'CType','bonferroni','Display','off');

% display stats
line([1 2],[34.5 34.5],'color',[0 0 0])
text(1.175, 36,'NS','color',[0 0 1])

line([2 3],[34 34],'color',[0 0 0])
text(2.175, 35.5,'NS','color',[1 0 0])

% ----------------------- save figure --------------------------
% saveas(gcf,[path_info.save_path,'decoding_figure2'])

%% more specific helper functions

function roi = get_roi_from_data(data)
% create binary roi mask from dFF data
    roi = data(:,:,1,1);
    roi(~isnan(roi)) = 1;
    roi(isnan(roi)) = 0;
end

function [reg_image, tform] = reg_image_by_roi(roi_ref,roi_fixed,data)
% image registration using binary ROI
    [optimizer, metric] = imregconfig('monomodal'); 
    tform = imregtform(roi_ref, roi_fixed,'similarity', optimizer, metric);

    data(isnan(data)) = 0; 
    reg_image = imwarp(data, tform, ...
        'OutputView', imref2d(size(roi_fixed)), 'FillValues', 0); 
end

function f_data = format_data(data, CL, CR, radius)
% format data for use with the ridge regression model
    f_data = zeros(size(data,4),size(CL,1),size(data,3));

    for i = 1:size(CL,1)
        regiony  = CL(i,2)-radius:CL(i,2)+radius;
        regionxL = CL(i,1)-radius:CL(i,1)+radius;
        regionxR = CR(i,1)-radius:CR(i,1)+radius;
        traces = squeeze(mean(mean(...
                    data(regiony,[regionxL,regionxR],:,:), ...
                    2),1));
        f_data(:,i,:) = permute(traces, [2 3 1]);
    end
end

function [subsampled_data, random_sample] = subsample_data(data1, data2, data3)
    [min_length, min_category] = min([size(data1,1) size(data2,1), size(data3,1)]);

    switch min_category
        case 1
        random_sample(1,:) = randi(size(data1,1),min_length,1);
        random_sample(2,:) = randi(size(data2,1),min_length,1);
        random_sample(3,:) = randi(size(data3,1),min_length,1);

        case 2
        random_sample(1,:) = randi(size(data1,1),min_length,1);
        random_sample(2,:) = randi(size(data2,1),min_length,1);
        random_sample(3,:) = randi(size(data3,1),min_length,1);

        case 3
        random_sample(1,:) = randi(size(data1,1),min_length,1);
        random_sample(2,:) = randi(size(data2,1),min_length,1);
        random_sample(3,:) = randi(size(data3,1),min_length,1);                       
    end

    subsampled_data{1} = data1(random_sample(1,:), :, :);
    subsampled_data{2} = data2(random_sample(2,:), :, :);
    subsampled_data{3} = data3(random_sample(3,:), :, :);
end

function effectsize = calc_effectsize(sample1,sample2)
    observeddifference = nanmean(sample1) - nanmean(sample2);
    effectsize = observeddifference / nanmean([std(sample1), std(sample2)]);
end