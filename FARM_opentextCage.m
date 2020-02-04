function expe = FARM_opentextCage(expe)

for i = 1:length(expe) % for each trial

    current_video_file = expe(i).filename; % current image name
    pos = findstr(current_video_file,'\'); % identify the location of "\"
    big_text_file = [current_video_file(1:pos(length(pos)-1)) 'TextFiles\headFix*.txt']; % set the expecte big text file for head-fixed data
    resu = dir(big_text_file); % check info about this text file
    big_text_file = [current_video_file(1:pos(length(pos)-1)) 'TextFiles\' resu.name]; % update the name of the text file
    current_video_file = current_video_file(max(pos)+2:end);
    disp(['checking ' current_video_file ' in ' big_text_file]) % display the identification of association between videofile and text

    try,
        evTable = readtable(big_text_file,'Delimiter','\t'); % read the text file
        evName = table2array(evTable(:,3)); % evName= name of the event

        no_file_found = 0; ii = 0;
        while no_file_found == 0 % while no association found, continue
            ii = ii + 1; % increment the event name
            texte = (char(evName(ii))); % convert in string
            q = findstr(texte,current_video_file); % find an event including the file name
            if length(q)==1 % if an event include file name (video.raw name)
                no_file_found = 1; % set that we found the location on text file
                no_stop_found = 0; jj = ii; % 
                debut = ii; % the starting is the current iteration
                while no_stop_found == 0 % now looking for the end
                    jj = jj + 1;
                    texte2 = (char(evName(jj)));
                    q2 = findstr(texte2,'complete');
                    q3 = findstr(texte2,'VideoEnd');
                    if length(q2)==1 || length(q3) == 1
                        no_stop_found = 1; 
                        fin = jj;
                    end
                end  
            end
        end

        expe(i).evTable = evTable(debut:fin,:);
    catch
        disp(i)
    end
end