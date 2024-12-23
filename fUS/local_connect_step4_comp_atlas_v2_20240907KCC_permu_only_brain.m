
%20241005 use only regions within brain

% load isocortex - cut off the region 

% plot overlap

Skcc=dir(fullfile('\\fenglab03\yiyun\20240820 fUS\','*kcc.mat'));

%S=cat(1,S1,S2);

%%
load('animalIDlist_final.mat');

hetorwt_kcc=cell(size(animalIDlistcom,1),3);



for si=1:numel(Skcc)
    
    filename=Skcc(si).name(2:end);
    if contains(filename,'_n')
        filename=filename(1:strfind(filename,'_n')-1);
    elseif ~isempty(strfind(filename,'r'))
        filename=filename(1:strfind(filename,'r')-1);
    end
    
    if ~isempty(strfind(filename,'k'))
        filename=filename(1:strfind(filename,'k')-1);
    end
    %print('filename')
    
    
    for ai=1:size(animalIDlistcom,1)
        animalID=animalIDlistcom{ai,1};
        if ~isempty(strfind(animalID,'l'))
            animalID=animalID(1:strfind(animalID,'l')-1);
        end
        if ~isempty(strfind(animalID,'r'))   
            animalID=animalID(1:strfind(animalID,'r')-1);            
        end
        animalID=strrep(animalID,'''','');
        
        %print('animalID')
        %animalID

        if strcmp(filename,animalID)==1
            hetorwt_kcc{si,1}=animalIDlistcom{ai,1};
            hetorwt_kcc{si,2}=animalIDlistcom{ai,2};
            hetorwt_kcc{si,3}=ai;

        end
        
        if strcmp(erase(filename,'_'),animalID)==1
            hetorwt_kcc{si,1}=animalIDlistcom{ai,1};
            hetorwt_kcc{si,2}=animalIDlistcom{ai,2};
            hetorwt_kcc{si,3}=ai;

        end
        
        if strcmp(erase(animalID,'_'),filename)==1
            hetorwt_kcc{si,1}=animalIDlistcom{ai,1};
            hetorwt_kcc{si,2}=animalIDlistcom{ai,2};
            hetorwt_kcc{si,3}=ai;

        end       
    end
    
end

%%
hetid_kcc=[];
wtid_kcc=[];

for si=1:size(hetorwt_kcc,1)
    
    if ~isempty(find(hetid==hetorwt_kcc{si,3}))
        hetid_kcc=[hetid_kcc;si];
    elseif ~isempty(find(wtid==hetorwt_kcc{si,3}))
        wtid_kcc=[wtid_kcc;si];
    end
            
end

%%
vqhet=[];
vqwt=[];
hetcount=0;
wtcount=0;
vqall=[];

for si=1:numel(Skcc)
    
%     if si<=35
%         load(strcat('\\fenglab03\yiyun\20230515 fUS no dex\local_connectivity_outliner_remove\',S(si).name),'vq_local');
%     else
%         load(strcat('\\fenglab03\yiyun\20240820 fUS\',S(si).name),'vq_local');
%     end

    load(strcat('\\fenglab03\yiyun\20240820 fUS\',Skcc(si).name),'KCC_local');
    
    if (sum(hetid_kcc==si)>0)

        if isempty(vqhet)
            vqhet=KCC_local;
        else
            vqhet=cat(4,vqhet,KCC_local);
        end
        hetcount=hetcount+1;

    elseif (sum(wtid_kcc==si)>0)
        if isempty(vqwt)
            vqwt=KCC_local;
        else
            vqwt=cat(4,vqwt,KCC_local);
        end
        wtcount=wtcount+1;

    end
    
    if (sum(hetid_kcc==si)>0)|(sum(wtid_kcc==si)>0)
        if isempty(vqall)
            vqall=KCC_local;
        else
            vqall=cat(4,vqall,KCC_local);
        end
    end
   
end

%%
vqwtav=mean(vqwt,4,'omitnan');
vqhetav=mean(vqhet,4,'omitnan');


%%
for i=1:5:41%21:30

    subplot(4,13,(i-1)/5+1)
    imagesc(permute(reshape(vqwtav(:,i,:),size(vqwtav,1),[]),[2,1]),[0.3 1]);
    %xlim([10 53]3
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('WT')
    end
    
    subplot(4,13,(i-1)/5+1+13)
    imagesc(permute(reshape(vqhetav(:,i,:),size(vqhetav,1),[]),[2,1]),[0.3 1]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('HET')
    end

    subplot(4,13,(i-1)/5+1+26)
    imagesc(permute(reshape((vqhetav(:,i,:)-vqwtav(:,i,:)),size(vqhetav,1),[]),[2,1]),[0 0.2]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('HET-WT')
    end

    subplot(4,13,(i-1)/5+1+39)
    imagesc(permute(reshape(-(vqhetav(:,i,:)-vqwtav(:,i,:)),size(vqhetav,1),[]),[2,1]),[0 0.2]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('WT-HET')
    end


end

%%

load('\\Fenglab03\Yiyun\20230410_fUS_seedmap\region_atlas_newly_sample_v5_241005resam_newreg.mat','vqTH','vqIsoCortex');
bdTH=bwperim(vqTH);
bdIsoCortex=bwperim(vqIsoCortex);

load('\\Fenglab03\Yiyun\20230410_fUS_seedmap\region_atlas_newly_sample_v5_241005MBresam_newreg.mat','vqMB');

%%
vqwtav2=vqwtav;

for i=1:size(vqwtav,2)
    for j=1:size(vqwtav,1)
        if ~isempty(find(vqIsoCortex(j,i,:)==1,1,'last'))
            vqwtav2(j,i,find(vqIsoCortex(j,i,:)==1,1,'last')+3:end)=0; %higher than ocrtex ==0
        else
            vqwtav2(j,i,find(vqMB(j,i,:)==1,1,'last')+3:end)=0; %higher than MB ==0
        end
        vqwtav2(j,i,1:10)=0;  %z lower than 10 ==0
    end
end

vqwtav2(:,[1:3,32:41],:)=0;

vqwtav2(68:81,:,:)=0;

%%
for i=1:1:41%21:30

    subplot(4,41,(i-1)/1+1)
    imagesc(permute(reshape(vqwtav2(:,i,:),size(vqwtav,1),[]),[2,1]),[0.5 0.9]);
    %xlim([10 53]3
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('WT')
    end
    
    subplot(4,41,(i-1)/1+1+41)
    imagesc(permute(reshape(vqhetav(:,i,:),size(vqhetav,1),[]),[2,1]),[0.5 0.9]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('HET')
    end

    subplot(4,41,(i-1)/1+1+41*2)
    imagesc(permute(reshape((vqwtav2(:,i,:)+0.3*vqTH(:,i,:)+0.3*vqIsoCortex(:,i,:)),size(vqhetav,1),[]),[2,1]),[0.3 0.8]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('HET-WT')
    end

    subplot(4,41,(i-1)/1+1+41*3)
    imagesc(permute(reshape((vqhetav(:,i,:)+0.3*vqTH(:,i,:)+0.3*vqIsoCortex(:,i,:)),size(vqhetav,1),[]),[2,1]),[0.3 0.8]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('WT-HET')
    end


end

%%
for i=1:3:41%21:30

    subplot(4,14,(i-1)/3+1)
    imagesc(permute(reshape(vqwtav2(:,i,:),size(vqwtav,1),[]),[2,1]),[0.5 0.9]);
    %xlim([10 53]3
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('WT')
    end
    
    subplot(4,14,(i-1)/3+1+14)
    imagesc(permute(reshape(vqhetav(:,i,:),size(vqhetav,1),[]),[2,1]),[0.5 0.9]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('HET')
    end

    subplot(4,14,(i-1)/3+1+14*2)
    imagesc(permute(reshape((vqwtav2(:,i,:)+0.3*vqTH(:,i,:)+0.3*vqIsoCortex(:,i,:)),size(vqhetav,1),[]),[2,1]),[0.3 0.8]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('HET-WT')
    end

    subplot(4,14,(i-1)/3+1+14*3)
    imagesc(permute(reshape((vqhetav(:,i,:)+0.3*vqTH(:,i,:)+0.3*vqIsoCortex(:,i,:)),size(vqhetav,1),[]),[2,1]),[0.3 0.8]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('WT-HET')
    end


end


%% remove non brain region in vqwt, vqhet and vqall


for k=1:size(vqwt,4)

    for i=1:size(vqwt,2)
        for j=1:size(vqwt,1)
            if ~isempty(find(vqIsoCortex(j,i,:)==1,1,'last'))
                vqwt(j,i,find(vqIsoCortex(j,i,:)==1,1,'last')+3:end,k)=0; %higher than ocrtex ==0
            else
                vqwt(j,i,find(vqMB(j,i,:)==1,1,'last')+3:end,k)=0; %higher than MB ==0
            end
            vqwt(j,i,1:10,k)=0;  %z lower than 10 ==0
        end
    end

    vqwt(:,[1:3,32:41],:,k)=0;

    vqwt(68:81,:,:,k)=0;

end

for k=1:size(vqhet,4)

    for i=1:size(vqhet,2)
        for j=1:size(vqhet,1)
            if ~isempty(find(vqIsoCortex(j,i,:)==1,1,'last'))
                vqhet(j,i,find(vqIsoCortex(j,i,:)==1,1,'last')+3:end,k)=0; %higher than ocrtex ==0
            else
                vqhet(j,i,find(vqMB(j,i,:)==1,1,'last')+3:end,k)=0; %higher than MB ==0
            end
            vqhet(j,i,1:10,k)=0;  %z lower than 10 ==0
        end
    end

    vqhet(:,[1:3,32:41],:,k)=0;

    vqhet(68:81,:,:,k)=0;

end

for k=1:size(vqall,4)

    for i=1:size(vqall,2)
        for j=1:size(vqall,1)
            if ~isempty(find(vqIsoCortex(j,i,:)==1,1,'last'))
                vqall(j,i,find(vqIsoCortex(j,i,:)==1,1,'last')+3:end,k)=0; %higher than ocrtex ==0
            else
                vqall(j,i,find(vqMB(j,i,:)==1,1,'last')+3:end,k)=0; %higher than MB ==0
            end
            vqall(j,i,1:10,k)=0;  %z lower than 10 ==0
        end
    end

    vqall(:,[1:3,32:41],:,k)=0;

    vqall(68:81,:,:,k)=0;

end

%%
%baseline

tstat_mat_base=zeros(size(vqwt,1),size(vqwt,2),size(vqwt,3));

for xi=1:size(vqwt,1)
    for yi=1:size(vqwt,2)
        for zi=1:size(vqwt,3)
            
            [h,p,c,st]=ttest2(vqwt(xi,yi,zi,:),vqhet(xi,yi,zi,:));
            tstat_mat_base(xi,yi,zi)=st.tstat;
            
        end
    end
end

save('KCC_ttest_baseline_brain_only.mat');

%%
clearvars -except vqall

%%

%permutation test

tstat_mat_per_max=[];

tic

for perti=1:5000
    
    hetwtlist=randperm(size(vqall,4));
    
    wtperi=hetwtlist(1:11);
    hetperi=hetwtlist(12:24);
    
    tstat_mat_per=zeros(size(vqall,1),size(vqall,2),size(vqall,3));
    for xi=1:size(vqall,1)
        for yi=1:size(vqall,2)
            for zi=1:size(vqall,3)

                [h,p,c,st]=ttest2(vqall(xi,yi,zi,wtperi),vqall(xi,yi,zi,hetperi));
                tstat_mat_per(xi,yi,zi)=st.tstat;

            end
        end
    end
    
    tstat_mat_per_max(perti)=max(reshape(tstat_mat_per,1,[]));
    
    if mod(perti,100)==1
        save('tstat_mat_per_max_mat_brain_only.mat')
    end
    
end


timeelapsed=toc

%%
%find voxel p <0.05

load('KCC_ttest_baseline_brain_only.mat');

load('tstat_mat_per_max_mat_brain_only.mat');


%%
sorted_tstat_mat_per_max=sort(tstat_mat_per_max,'descend');


%find voxel p<0.05, 0.05*5000=250 250+1
[a,b,c]=ind2sub(size(tstat_mat_base),find(tstat_mat_base>=sorted_tstat_mat_per_max(251)));

%%
%%
load('\\fenglab03\yiyun\20230410_fUS_seedmap\region_atlas_newly_sample_v5_231011resam.mat','vqACC','vqPL','vqIL','vqMD','vqSTR','vqRSP','vqTRN','vqRE','vqSTRd','vqSTRv');
%load('\\fenglab03\yiyun\20230410_fUS_seedmap\region_atlas_newly_sample_v5re_strd_strv.mat',);
load('\\Fenglab03\Yiyun\20230410_fUS_seedmap\region_atlas_newly_sample_v5_240823resam_newreg.mat','vqSSp','vqDORsm','vqDORpm','vqHPF','vqMO','vqHY','vqSSs','vqPTLp');

vqDORpml=vqDORpm;
vqDORpmr=vqDORpm;

vqMDr=vqMD;
vqMDr(:,1:20,:)=0;

vqDORpml(:,1:20,:)=0;
vqDORpmr(:,21:41,:)=0;

bdMD=bwperim(vqMD);
bdTRN=bwperim(vqTRN);
bdACC=bwperim(vqACC);
bdPL=bwperim(vqPL);
bdIL=bwperim(vqIL);
bdRSP=bwperim(vqRSP);
%bdS1=bwperim(vqS1);
bdSTR=bwperim(vqSTR);

bdSTRd=bwperim(vqSTRd);
bdSTRv=bwperim(vqSTRv);
bdRE=bwperim(vqRE);


bdMO=bwperim(vqMO);
bdHPF=bwperim(vqHPF);
bdHY=bwperim(vqHY);
bdDORsm=bwperim(vqDORsm);
bdDORpm=bwperim(vqDORpm);

bdDORpml=bwperim(vqDORpml);
bdDORpmr=bwperim(vqDORpmr);


for i=3:3:41%21:30

    subplot(4,13,(i-3)/3+1)
    imagesc(permute(reshape(tstat_mat_base(:,i,:),size(tstat_mat_base,1),[]),[2,1]),[-2.5 4.5]);
    xlim([0 67])
    ylim([10 74])
    set(gca, 'YDir','normal');
    if i==3
    title('t-map')
    end
    
    subplot(4,13,(i-3)/3+13+1)
    imagesc(permute(reshape((tstat_mat_base(:,i,:)>=sorted_tstat_mat_per_max(251)),size(tstat_mat_base,1),[]),[2,1]));
    xlim([0 67])
    ylim([10 74])
    set(gca, 'YDir','normal');
    
    title(strcat(num2str(i),'th slice'))
   
    if i==21
    title('midline')
    end

    subplot(4,13,(i-3)/3+13*2+1)
    imagesc(permute(reshape(tstat_mat_base(:,i,:)+5*(~bdMD(:,i,:)).*(~bdTRN(:,i,:)).*(~bdACC(:,i,:)).*(~bdPL(:,i,:)).*(~bdIL(:,i,:)).*(~bdRSP(:,i,:)).*(~bdSTRd(:,i,:)).*(~bdSTRv(:,i,:)).*(~bdRE(:,i,:)).*(~bdMO(:,i,:)).*(~bdHPF(:,i,:)).*(~bdHY(:,i,:)),size(tstat_mat_base,1),[]),[2,1]),[2 8]);
    xlim([0 67])
    ylim([10 74])
    set(gca, 'YDir','normal');
    if i==3
    title('t-map')
    end

end

set(gcf,'renderer','painters')
%%
for i=15

    subplot(3,1,1)
    imagesc(permute(reshape(tstat_mat_base(:,i,:),size(tstat_mat_base,1),[]),[2,1]),[-5 7]);
    %xlim([10 53]3
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==1
    title('t-map')
    end
    
    subplot(3,1,2)
    imagesc(permute(reshape((tstat_mat_base(:,i,:)>=sorted_tstat_mat_per_max(251)),size(tstat_mat_base,1),[]),[2,1]));
    %xlim([10 53]3
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==1
    title('t-map')
    end

    subplot(3,1,3)
    imagesc(permute(reshape((vqDORpmr(:,i,:)),size(tstat_mat_base,1),[]),[2,1]),[0 1]);
    %xlim([10 53]3
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==1
    title('t-map')
    end

end