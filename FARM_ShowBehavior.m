function FARM_ShowBehavior(expe)
try
    
clear SD MM G RMS
for i = 1:length(expe)
    for j = 1:length(expe(i).be)
        SD(:,:,i,j) = expe(i).be(j).SD;
        MM(:,:,i,j) = expe(i).be(j).MM;
        G(:,:,i,j) = expe(i).be(j).G;
        RMS(:,:,i,j) = expe(i).be(j).RMS;
    end
end

SDt = [];
Gt = [];
MMt = [];
RMSt = [];
for i = 1:size(SD,4) % for each cam    
    SDi = SD(:,:,:,i); % list of SD
    SDi = SDi(:,:,squeeze(sum(sum(SDi,1),2))>0); % keep only those with data, not zero
    SDt = cat(2,SDt,OIA_n(mean(SDi,3))); % average, normalize from 0 to 1 and concactenate.
    
    MMi = MM(:,:,:,i); % same for MM
    MMi = MMi(:,:,squeeze(sum(sum(MMi,1),2))>0); 
    MMt = cat(2,MMt,OIA_n(mean(MMi,3))); 

    Gi = G(:,:,:,i); % same for MM
    Gi = Gi(:,:,squeeze(sum(sum(Gi,1),2))>0); 
    Gt = cat(2,Gt,OIA_n(mean(Gi,3))); 
    
    RMSi = RMS(:,:,:,i); % same for MM
    RMSi = RMSi(:,:,squeeze(sum(sum(RMSi,1),2))>0); 
    RMSt = cat(3,RMSt,OIA_n(mean(RMSi,3))); 
    
end

RESU = cat(1,MMt,SDt,Gt,RMSt);
imshow(RESU)

end