function [expe,varOpen] = FARM_Process_beh(expe,varOpen)

try, disp(['              varOpen.varProc.gsr=' num2str(varOpen.varProc.gsr)]); 
    catch, varOpen.varProc.gsr = 0; disp('(not defined) varOpen.varProc.gsr=0'); end
try, disp(['              varOpen.varProc.f=' num2str(varOpen.varProc.f)]); 
    catch, varOpen.varProc.f = []; disp('(not defined) varOpen.varProc.f=[]'); end
try, disp(['              varOpen.varProc.Ca_flipLR=' num2str(varOpen.varProc.Ca_flipLR)]); 
    catch, varOpen.varProc.Ca_flipLR = 0; disp('(not defined) varOpen.varProc.Ca_flipLR=0'); end
try, disp(['              varOpen.varProc.BinningCa=' num2str(varOpen.varProc.BinningCa)]); 
    catch, varOpen.varProc.BinningCa = 1; disp('(not defined) varOpen.varProc.BinningCa=1'); end
try, disp(['              varOpen.varProc.tb=' num2str(varOpen.varProc.tb)]); 
    catch, varOpen.varProc.tb = 1; disp('(not defined) varOpen.varProc.tb=1'); end
try, disp(['              varOpen.varProc.blueCorr=' num2str(varOpen.varProc.blueCorr)]); 
    catch, varOpen.varProc.blueCorr = 1; disp('(not defined) varOpen.varProc.blueCorr=1'); end
try, disp(['              varOpen.varProc.run_list=' num2str(varOpen.varProc.run_list)]); 
    catch, varOpen.varProc.run_list = []; disp('(not defined) varOpen.varProc.run_list=[]'); end
try, disp(['              varOpen.varProc.useBe1=' num2str(varOpen.varProc.useBe1)]); 
    catch, varOpen.varProc.useBe1 = 0; disp('(not defined) varOpen.varProc.useBe1=0'); end
try, disp(['              varOpen.varProc.useBe2=' num2str(varOpen.varProc.useBe2)]); 
    catch, varOpen.varProc.useBe2 = 0; disp('(not defined) varOpen.varProc.useBe2=0'); end
try, disp(['              varOpen.varProc.useBe3=' num2str(varOpen.varProc.useBe3)]); 
    catch, varOpen.varProc.useBe3 = 0; disp('(not defined) varOpen.varProc.useBe3=0'); end
                                        
% 2. Identify the trial who will be used ----------------------------------

roi = OIA_n(expe(1).roi); % roi = roi of the first run
bx = expe(1).bx; % bregma coordinates
by = expe(1).by;
if varOpen.varProc.BinningCa ~=1, % if binning calcium
    roi = imresize(roi,varOpen.varProc.BinningCa,'nearest'); 
    bx = round(bx*varOpen.varProc.BinningCa);
    by = round(by*varOpen.varProc.BinningCa);
end
if varOpen.varProc.Ca_flipLR == 1, 
    roi = fliplr(roi); 
    by = size(roi,1)-by;
end;
varOpen.varOutput.bx = bx; % transfer bregma
varOpen.varOutput.by = by; % transfer bregma


roinan = roi;
try, roinan(roinan==0) = NaN;
end
varOpen.varOutput.roi = roi;
varOpen.varOutput.roinan = roinan;

sBe = varOpen.varProc.useBe1 + varOpen.varProc.useBe2 + varOpen.varProc.useBe3;
if varOpen.varProc.useBe1 == 0, varOpen.varProc.useBe1 == NaN; end;
if varOpen.varProc.useBe2 == 0, varOpen.varProc.useBe2 == NaN; end;
if varOpen.varProc.useBe3 == 0, varOpen.varProc.useBe3 == NaN; end;

is_there_behavior = zeros(3,length(expe));
for i = 1:3 % check for each possible behavior view
    for j = 1:length(expe) % check for each run
        try, 
            is_there_behavior(i,j) = size(expe(j).be(i).J,3); % size of behavior file, if any
        end
    end
end
is_there_behavior = is_there_behavior>1; % only consider a good behavior file if more than one frame
Be1 = varOpen.varProc.useBe1 .* is_there_behavior(1,:); % determine which beh cameras to take into account 
Be2 = varOpen.varProc.useBe2 .* is_there_behavior(2,:);
Be3 = varOpen.varProc.useBe3 .* is_there_behavior(3,:);
if sBe == 0, Be = ones(1,length(expe)); % if no behavior requested, use all trial
else, Be = (nansum([Be1; Be2; Be3])/sBe)>=1; end; % else check only trial where all the condition are met
varOpen.varOutput.sBe = sBe; % report number of behavior used

if length(varOpen.varProc.run_list) ==0, % if no list manually established, take all trial who met the condition
    varOpen.varProc.run_list = find(Be==1); 
    run_list = varOpen.varProc.run_list;
else, % else, only considere those who met the condition
    run_list = varOpen.varProc.run_list;
    run_test = zeros(1,length(expe));
    for i = 1:length(run_list)
        run_test(run_list(i)) = 1;
    end
    Be2 = run_test .* Be;
    varOpen.varProc.run_list = find(Be2==1); 
end; 


% 3. Concatenation procedure ----------------------------------------------

REF = []; % reference image (fluorescence)
REF2 = []; % reference image (reflectance)
II = []; % concatenated calcium matrix (fluo)
JJ = []; % concatenated behavior matrix (if any)
SD = []; G = []; RMS = []; MM = [];
jonction = [];
SFr = [];
for i = 1:length(varOpen.varProc.run_list) % for each trial considered
    disp(['---- processing run ' num2str(varOpen.varProc.run_list(i)) ' -------------'])
    
    
    evTable = char(table2array(expe(varOpen.varProc.run_list(i)).evTable(:,3))); % Text file description
    evTime = table2array(expe(varOpen.varProc.run_list(i)).evTable(:,2)); % Time corresponding to every text file event 
    is_LED_tag = 0 ;
    for k = 1:size(evTable,1) % go through all entries to find BrainLEDON and OFF times 
            q = findstr(evTable(k,:),'BrainLEDON'); % determine BrainLEDON and OFF times for every expe from text file
            if length(q)==1
                BrainLEDON = evTime(k)-evTime(1); % subtracting actual time by the first time in the table 
                is_LED_tag = is_LED_tag+1; 
            else
            end
            f = findstr(evTable(k,:),'BrainLEDOFF'); 
            if length(f)==1
                BrainLEDOFF = evTime(k)-evTime(1);
                is_LED_tag = is_LED_tag+1; 
            else
            end
    end
    
    BrainLEDSUB = BrainLEDOFF - BrainLEDON; % Total length of time being taken into account 
    
    if is_LED_tag == 2 % only if BrainLEDON (+1) and BrainLEDOFF (+1) are present 
            realSF = size(expe(varOpen.varProc.run_list(i)).I,3)/(BrainLEDSUB);
            expe(varOpen.varProc.run_list(i)).realSF = realSF; % real sampling frequency of calcium data for each expe 
            realSF_list(i) = realSF; 
    else
    end 
    
    varOpen.varProc.events(varOpen.varProc.run_list(i)).BrainLEDOFF = BrainLEDOFF; % save individual expe data in varOpen 
    varOpen.varProc.events(varOpen.varProc.run_list(i)).BrainLEDON = BrainLEDON; 
    varOpen.varProc.events(varOpen.varProc.run_list(i)).BrainLEDSUB = BrainLEDSUB;
        
    If = single(expe(varOpen.varProc.run_list(i)).I); % pick up calcium data from every expe 
    expe(varOpen.varProc.run_list(i)).REF = []; % clear the reference image
    if varOpen.varProc.BinningCa ~=1, If = imresize(If,varOpen.varProc.BinningCa,'nearest'); end; % calcium binning
    if varOpen.varProc.Ca_flipLR == 1, If = fliplr(If); end; % if flipping

    for j = 1:size(If,4) % for each channel (=color)
        expe(varOpen.varProc.run_list(i)).REF(:,:,j) = mean(If(:,:,:,j),3); % Reference image (= time average)
        if varOpen.varProc.gsr == 1, % if GSR, do it for each frame
            clear v, v.method = 1; v.dc = 0; v.roi = roi; v.disp = 0;
            If(:,:,:,j) = OIA_SA_globalactivity(If(:,:,:,j),v); 
        else % if no GSR, just substract DC
            If(:,:,:,j) = OIA_addDC(If(:,:,:,j),-expe(varOpen.varProc.run_list(i)).REF(:,:,j));
        end
        if length(varOpen.varProc.f) ~= 0 % if filtering
            clear tf, tf.method = 3; tf.lp = varOpen.varProc.f(1); tf.hp = varOpen.varProc.f(2); 
            tf.SF = SF; tf.disp = 0;
            If(:,:,:,j) = OIA_SA_temporalfilter(If(:,:,:,j),tf);
        end
        current_roinan = imresize(roinan,size(expe(varOpen.varProc.run_list(i)).REF(:,:,j)));
        If(:,:,:,j) = 100*OIA_multiDC(If(:,:,:,j),1./expe(varOpen.varProc.run_list(i)).REF(:,:,j)); % calculate DF/F
        If(:,:,:,j) = OIA_multiDC(If(:,:,:,j),current_roinan);
        %expe(varOpen.varProc.run_list(i)).If(:,:,:,j) = OIA_multiDC(If(:,:,:,j),current_roinan); % store DF/F with mask
    end
   
    REF = cat(3,REF,single(expe(varOpen.varProc.run_list(i)).REF(:,:,1))); % concatenate Reference (fluo)
    try, REF2 = cat(3,REF2,single(expe(varOpen.varProc.run_list(i)).REF(:,:,2))); end % concatenate Reference (reflectance, if any)
    
    JJij = [];
    SDij = []; MMij = []; RMSij = []; Gij = []; 
    % include behavior if any (and give feedback) - cam 1
    if varOpen.varProc.useBe1 == 1 
        disp(' > adding behavior cam 1')
        JJi = expe(varOpen.varProc.run_list(i)).be(1).J; % pick up the behavior run; 
        SDij = cat(2,SDij,expe(varOpen.varProc.run_list(i)).be(1).SD); % standard dev
        Gij = cat(2,Gij,expe(varOpen.varProc.run_list(i)).be(1).G); % avg gradient 
        MMij = cat(2,MMij,expe(varOpen.varProc.run_list(i)).be(1).MM); % avg over time 
        RMSij = cat(2,RMSij,expe(varOpen.varProc.run_list(i)).be(1).RMS); % root mean square
        realSF_cam = size(expe(varOpen.varProc.run_list(i)).be(1).J,3)/varOpen.varProc.events(varOpen.varProc.run_list(i)).BrainLEDSUB; % calculating sampling frequency for the behavioural data for every expe  

        if varOpen.varProc.tb ~= 1
            JJi = OIA_tempobin(JJi,1/varOpen.varProc.tb); % temporal binning
        end
        JJij = cat(2,JJij,JJi); % concatenate behavioural videos to have them display side to side 
        
    end
    
    % include behavior if any (and give feedback) - cam 2
    if varOpen.varProc.useBe2 == 1
        disp(' > adding behavior cam 2')
        JJi = expe(varOpen.varProc.run_list(i)).be(2).J; % pick up the behavior run
        SDij = cat(2,SDij,expe(varOpen.varProc.run_list(i)).be(2).SD);
        Gij = cat(2,Gij,expe(varOpen.varProc.run_list(i)).be(2).G);
        MMij = cat(2,MMij,expe(varOpen.varProc.run_list(i)).be(2).MM);
        RMSij = cat(2,RMSij,expe(varOpen.varProc.run_list(i)).be(2).RMS);
        if sum(size(JJij)) ~= 0 % if data already exists 
            if size(JJi,3) ~= size(JJij,3) % resizing it in the time dim so data can be concatenated 
                JJi1= reshape(JJi,[size(JJi,1)*size(JJi,2), size(JJi,3)]); 
                JJi2= imresize(JJi1,[size(JJi1,1), size(JJij,3)]); 
                JJi= reshape(JJi2,[size(JJi,1), size(JJi,2), size(JJi2,2)]); 
            end 
        else 
            realSF_cam= size(expe(varOpen.varProc.run_list(i)).be(2).J,3)/varOpen.varProc.events(varOpen.varProc.run_list(i)).BrainLEDSUB; 
            if varOpen.varProc.tb ~= 1
            JJi = OIA_tempobin(JJi,1/varOpen.varProc.tb); % temporal binning
            end
        end 
 
        JJij = cat(2,JJij,JJi); % concaternate in X-axis the behavior view
    end
    
    % include behavior if any (and give feedback) - cam 3
    if varOpen.varProc.useBe3 == 1
        disp(' > adding behavior cam 3')
        JJi = expe(varOpen.varProc.run_list(i)).be(3).J; % pick up the behavior run
        SDij = cat(2,SDij,expe(varOpen.varProc.run_list(i)).be(3).SD);
        Gij = cat(2,Gij,expe(varOpen.varProc.run_list(i)).be(3).G);
        MMij = cat(2,MMij,expe(varOpen.varProc.run_list(i)).be(3).MM);
        RMSij = cat(2,RMSij,expe(varOpen.varProc.run_list(i)).be(3).RMS);
        if (varOpen.varProc.useBe2 + varOpen.varProc.useBe1) == 2 
            if size(JJi,3) ~= size(JJij,3) % resizing it in the time dim so data can be concatenated 
                JJi1= reshape(JJi,[size(JJi,1)*size(JJi,2), size(JJi,3)]); 
                JJi2= imresize(JJi1,[size(JJi1,1), size(JJij,3)]); 
                JJi= reshape(JJi2,[size(JJi,1), size(JJi,2), size(JJi2,2)]); 
            end 
        elseif (varOpen.varProc.useBe2 + varOpen.varProc.useBe1) == 1
            if size(JJi,3) ~= size(JJij,3)
                JJi1= reshape(JJi,[size(JJi,1)*size(JJi,2), size(JJi,3)]); 
                JJi2= imresize(JJi1,[size(JJi1,1), size(JJij,3)]); 
                JJi= reshape(JJi2,[size(JJi,1), size(JJi,2), size(JJi2,2)]); 
            end 
        else % if nothing in concatenated data
            if varOpen.varProc.tb ~= 1
            JJi = OIA_tempobin(JJi,1/varOpen.varProc.tb); % temporal binning
            end
            realSF_cam = size(expe(varOpen.varProc.run_list(i)).be(3).J,3)/varOpen.varProc.events(varOpen.varProc.run_list(i)).BrainLEDSUB; 
        end 
%         if size(JJi,3)~=size(If,3) % if not the same length than calcium... 
%             JJi = OIA_tempobin(JJi,size(If,3)/size(JJi,3)); % ...adjust the number of frame to fit
%         end
        JJij = cat(2,JJij,JJi); % concaternate in X-axis the behavior view
    end
        
    expe(varOpen.varProc.run_list(i)).realSF_cam= realSF_cam; % real sampling frequency of behavioural 
    expe(varOpen.varProc.run_list(i)).Beh_cam = JJij; % save expe behavioural videos to their own expe field  
    SD = cat(3,SD,SDij);
    G = cat(3,G,Gij);
    RMS = cat(3,RMS,RMSij);
    MM = cat(3,MM,MMij);
    
    JJ = cat(3,JJ,JJij); % concaternate in T-axis the behavior view (of all cam)
    If = OIA_tempobin(If,1/varOpen.varProc.tb); % If: calcium data for each expe 
    II = cat(3,II,If); % concatenation of fluo (and reflectance)
    jonction(run_list(i)) = size(If,3); % size of the concatenated data
    
    expe(varOpen.varProc.run_list(i)).If = If; % individual calcium videos 
    expe(varOpen.varProc.run_list(i)).REFbe = JJij(:,:,1); 
    
    SFr(i) = expe(varOpen.varProc.run_list(i)).SF .* size(If,3)/size(expe(varOpen.varProc.run_list(i)).I,3);
    
    
    
end
varOpen.varOutput.REF = mean(REF,3); % average reference map (basal fluo)
try, varOpen.varOutput.REF2 = mean(REF2,3); end % average reference map (basal reflectance)


varOpen.varOutput.I = II; % transfering concatenated fluo (and refl)
varOpen.varOutput.Be = JJ; % transfering concatenatedbehavior
varOpen.varOutput.jonction = jonction;
varOpen.varOutput.SD = mean(SD,3);
varOpen.varOutput.G = mean(G,3);
varOpen.varOutput.RMS = mean(RMS,3);
varOpen.varOutput.MM = mean(MM,3);
varOpen.varOutput.SFr = mean(SFr);
% 4. check the text data --------------------------------------------------
LEDOFF = 0; % concatenated offset between runs
event_list = []; % time of licks 
event_list_reward = []; % time of reward
event_list_task = []; % time of cue given 
event_list_task_info = []; % description of type of cue 
realSF_list = []; % real sampling frequency of calcium 
time_event_expe= []; 
event_reward= [];
time_task_expe= []; 

disp('-------------------------------------------------------------------')
for i = 1:length(varOpen.varProc.run_list) % for each run
    disp([' > formating text information (trial ' num2str(varOpen.varProc.run_list(i)) ')'])
    try,
        range= [] ; 
        evTable = expe(varOpen.varProc.run_list(i)).evTable; % pick-up table
        evName = table2array(evTable(:,3)); % event name row
        evTime = table2array(evTable(:,2)); % real time row
        if iscell(evTime) == 1, evTime = str2num(cell2mat(evTime)); end

        time_event = []; % lick 
        time_reward = []; % reward
        time_task = []; % cue 
        time_task_info = []; % type of cue 
 
        is_LED_tag = 0; % reporter to check if LED on Tag

        for j = 1:length(evName) % for each event name of the current file
            texte = (char(evName(j))); % check the string for:
            q = findstr(texte,'lick:'); % ... lick
            if length(q)==1
                time_event = [time_event evTime(j)-evTime(1)]; % report this time relative to the initial time
            end

            q = findstr(texte,'reward'); 
            if length(q)==1
                time_reward = [time_reward evTime(j)-evTime(1)];
            end

            q = findstr(texte,'lickWitholdTime'); if length(q)==1 % cue 
                time_task = [time_task evTime(j)-evTime(1)];
                time_task_info_i = evName(j);
                time_task_info = [time_task_info time_task_info_i];

            end

            q = findstr(texte,'BrainLEDON'); if length(q)==1
                BrainLEDON = evTime(j)-evTime(1); 
                is_LED_tag = is_LED_tag+1; end
            q = findstr(texte,'BrainLEDOFF'); if length(q)==1
                BrainLEDOFF = evTime(j)-evTime(1);
                is_LED_tag = is_LED_tag+1; 
            end

        end

        
        if is_LED_tag==2 % only if BrainLEDON (+1) and BrainLEDOFF (+1) are present
            realSF = size(expe(varOpen.varProc.run_list(i)).I,3)/(BrainLEDOFF-BrainLEDON); % calcium SF 
            expe(varOpen.varProc.run_list(i)).realSF = realSF;
            realSF_list(i) = realSF;

            time_event = time_event - BrainLEDON; % all the event time are substracted to be relative the LED ON period
            time_event = time_event(time_event>0); % exclude data before BrainLEDON time 
            time_reward = time_reward - BrainLEDON; % same for reward...
            time_reward = time_reward(time_reward>0); 
            time_task = time_task - BrainLEDON; % same for task...
            time_task = time_task(time_task>0); 

            LEDOFFi = BrainLEDOFF-BrainLEDON; % time accounted for 

            time_event = time_event(time_event<LEDOFFi); % keep only event before LED-OFF
            time_reward = time_reward(time_reward<LEDOFFi); 
            time_task = time_task(time_task<LEDOFFi);
            
            range = [range LEDOFF] ; % indicate which times each expe is responsible for 
            event_list = [event_list time_event+LEDOFF]; % contactenate event tag lick
            time_event_expe= [ time_event_expe repelem(varOpen.varProc.run_list(i), length(time_event))];
            event_list_reward = [event_list_reward time_reward+LEDOFF]; % contactenate event tag lick
            event_reward = [ event_reward repelem(varOpen.varProc.run_list(i), length(time_reward))];
            event_list_task = [event_list_task time_task+LEDOFF]; % contactenate event tag lick
            time_task_expe = [time_task_expe repelem(varOpen.varProc.run_list(i),length(time_task))] ;
            event_list_task_info = [event_list_task_info time_task_info];
            varOpen.varProc.events(varOpen.varProc.run_list(i)).event_lick = time_event; % time of lick in individual expe 
            varOpen.varProc.events(varOpen.varProc.run_list(i)).lick_LED = time_event+LEDOFF; % take into account the concatenated times of the trials for licks 
            varOpen.varProc.events(varOpen.varProc.run_list(i)).event_reward = time_reward; % time of reward given in individual expe
            varOpen.varProc.events(varOpen.varProc.run_list(i)).reward_LED = time_reward+LEDOFF; % concatenated time
            varOpen.varProc.events(varOpen.varProc.run_list(i)).event_task = time_task; % time that GO trial cue is given in every expe 
            varOpen.varProc.events(varOpen.varProc.run_list(i)).task_LED = time_task +LEDOFF; % concatenated time
            varOpen.varProc.events(varOpen.varProc.run_list(i)).event_task_info= time_task_info; 

            LEDOFF = LEDOFF + LEDOFFi; % LED OFF 
            range = [range LEDOFF] ; 
            varOpen.varProc.events(varOpen.varProc.run_list(i)).time_window = range; 

        else
            disp('No BrainLEDON and OFF tags')
        end
         
    catch
        disp('no behavior text file for that run')
    end

end

varOpen.varOutput.LEDOFF = LEDOFF; % concatenated offset between runs
varOpen.varOutput.event_list = event_list; % concatenated lick times 
varOpen.varOutput.event_list_expe = time_event_expe; 
varOpen.varOutput.event_list_reward = event_list_reward; % reward
varOpen.varOutput.event_list_reward_expe = event_reward;
varOpen.varOutput.event_list_task = event_list_task; % go trial tasks 
varOpen.varOutput.event_list_task_expe = time_task_expe;
varOpen.varOutput.event_list_task_info = event_list_task_info;
varOpen.varOutput.realSF_list = realSF_list;

% 7. This part generate a video of the trial/trial variability of the view for calcium but also for the behavior

figure
REFv = []; % this will be the matrix used to feed OIA_new_video function
for i = 1:length(varOpen.varProc.run_list) % for each trial
    REFca = OIA_n(expe(varOpen.varProc.run_list(i)).REF(:,:,1)); % first, REFca is the reference image of calcium
    if size(expe(varOpen.varProc.run_list(i)).REF,3) == 2 % if there is a second channel = reflectance
        REFca = cat(1,REFca,OIA_n(expe(varOpen.varProc.run_list(i)).REF(:,:,2)));
    end

    try,
        REFbe = OIA_n(single(expe(varOpen.varProc.run_list(i)).REFbe)); % try to load reference image of behavior, camera 1
        REFca = imresize(REFca,size(REFbe,1)/size(REFca,1),'nearest'); % change de size of calcium to fit with behavior size
        REFv = cat(3,REFv,cat(2,REFca,REFbe)); % put side by side, calcium and behavior
    catch
       REFv = cat(3,REFv,REFca); % if no behavior, just add calcium
    end 
end

varOpen.varOutput.REFv = REFv; % transfer diplay matrix
varOpen.varOutput.L = expe(1).L; % transfer width
