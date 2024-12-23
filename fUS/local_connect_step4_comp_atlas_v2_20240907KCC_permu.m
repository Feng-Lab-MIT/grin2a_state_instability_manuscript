%use this after running the step1
%230911 add RE,vSTR, dSTR
%for permutation test

%S1=dir(fullfile('\\fenglab03\yiyun\20230515 fUS no dex\local_connectivity_outliner_remove\','*nodexlocalcon_gf1011.mat'));
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

save('KCC_ttest_baseline.mat');

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
        save('tstat_mat_per_max_mat.mat')
    end
    
end


timeelapsed=toc

%%

%find voxel p <0.05

sorted_tstat_mat_per_max=sort(tstat_mat_per_max,'descend');


%find voxel p<0.05, 0.05*5000=250 250+1
[a,b,c]=ind2sub(size(tstat_mat_base),find(tstat_mat_base>=sorted_tstat_mat_per_max(251)));

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


for i=3:4:41%21:30

    subplot(4,10,(i-3)/4+1)
    imagesc(permute(reshape(tstat_mat_base(:,i,:),size(tstat_mat_base,1),[]),[2,1]),[-3 5]);
    %xlim([10 53]3
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==1
    title('t-map')
    end
    
    subplot(4,10,(i-3)/4+10+1)
    imagesc(permute(reshape((tstat_mat_base(:,i,:)>=sorted_tstat_mat_per_max(251)),size(tstat_mat_base,1),[]),[2,1]));
    %xlim([10 53]3
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==1
    title('t-map')
    end

    subplot(4,10,(i-3)/4+10*2+1)
    imagesc(permute(reshape((~bdMD(:,i,:)).*(~bdTRN(:,i,:)).*(~bdACC(:,i,:)).*(~bdPL(:,i,:)).*(~bdIL(:,i,:)).*(~bdRSP(:,i,:)).*(~bdSTRd(:,i,:)).*(~bdSTRv(:,i,:)).*(~bdRE(:,i,:)).*(~bdMO(:,i,:)).*(~bdHPF(:,i,:)).*(~bdHY(:,i,:)),size(tstat_mat_base,1),[]),[2,1]),[0 1]);
    %xlim([10 53]3
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==1
    title('t-map')
    end

end

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