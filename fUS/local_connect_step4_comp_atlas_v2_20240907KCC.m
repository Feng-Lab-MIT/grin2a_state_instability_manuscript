%use this after running the step1
%230911 add RE,vSTR, dSTR

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
vq_all=[];

vqhet=[];
vqwt=[];
hetcount=0;
wtcount=0;

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
   
end


%%
%get only the brain region
load('\\Fenglab03\Yiyun\20230410_fUS_seedmap\region_atlas_newly_sample_v5_241005resam_newreg.mat','vqTH','vqIsoCortex'); %cortex
load('\\Fenglab03\Yiyun\20230410_fUS_seedmap\region_atlas_newly_sample_v5_241005MBresam_newreg.mat','vqMB'); %olfactory lobe

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

%%
vqwtav=mean(vqwt,4,'omitnan');
vqhetav=mean(vqhet,4,'omitnan');


%%
for i=1:3:41%21:30

    subplot(4,14,(i-1)/3+1)
    imagesc(permute(reshape(vqwtav(:,i,:),size(vqwtav,1),[]),[2,1]),[0.5 0.8]);
    xlim([0 67])
    ylim([10 74])
    set(gca, 'YDir','normal');
    if i==11
    title('WT')
    end
    
    subplot(4,14,(i-1)/3+1+14)
    imagesc(permute(reshape(vqhetav(:,i,:),size(vqhetav,1),[]),[2,1]),[0.5 0.8]);
    xlim([0 67])
    ylim([10 74])
    set(gca, 'YDir','normal');
    if i==11
    title('HET')
    end

    subplot(4,14,(i-1)/3+1+28)
    imagesc(permute(reshape((vqhetav(:,i,:)-vqwtav(:,i,:)),size(vqhetav,1),[]),[2,1]),[0 0.15]);
    xlim([0 67])
    ylim([10 74])
    set(gca, 'YDir','normal');
    if i==11
    title('HET-WT')
    end

    subplot(4,14,(i-1)/3+1+42)
    imagesc(permute(reshape(-(vqhetav(:,i,:)-vqwtav(:,i,:)),size(vqhetav,1),[]),[2,1]),[0 0.15]);
    xlim([0 67])
    ylim([10 74])
    set(gca, 'YDir','normal');
    if i==11
    title('WT-HET')
    end


end

%%

load('\\fenglab03\yiyun\20230410_fUS_seedmap\region_atlas_newly_sample_v5_231011resam.mat','vqACC','vqPL','vqIL','vqMD','vqSTR','vqRSP','vqTRN','vqRE','vqSTRd','vqSTRv');
%load('\\fenglab03\yiyun\20230410_fUS_seedmap\region_atlas_newly_sample_v5re_strd_strv.mat',);
load('\\Fenglab03\Yiyun\20230410_fUS_seedmap\region_atlas_newly_sample_v5_240823resam_newreg.mat','vqSSp','vqDORsm','vqDORpm','vqHPF','vqMO','vqHY','vqSSs','vqPTLp');


%vqACC=vqACC;
%vqACCa(1:40,:,:)=0;
%vqACCa(45:65,:,:)=0;%remove where brain disappear
%vqPL(45:65,:,:)=0;%remove where brain disappear

% vqACCal=vqACCa;
%vqACCar=vqACCa;
% 
% vqACCal(:,1:30,:)=0;
vqACCr=vqACC;
vqACCr(:,21:41,:)=0;

vqILr=vqIL;
vqILr(:,21:41,:)=0;

vqPLr=vqPL;
vqPLr(:,21:41,:)=0;

% 1. motor cortex (MO)
% 2. RSP retrosplenial area
% 3. dorsal striatum dSTR
% 4. ventral striatum vSTR
% 5. hippocampus HFP
% 6. Hypothalamus HY
% 7. DORsm Thalamus, sensory-motor cortex related
% 8. DORpm Thalamus, polymodal association cortex

% PTLp  Posterior parietal association area
%  midline thalamus 
%  lateral thalamus
%  ventral thalamus 
%  RE

vqMOl=vqMO;
vqMOr=vqMO;

vqMOl(:,1:20,:)=0;
vqMOr(:,21:41,:)=0;

vqSTRvl=vqSTRv;
vqSTRvr=vqSTRv;

vqSTRvl(:,1:20,:)=0;
vqSTRvr(:,21:41,:)=0;

vqSTRdl=vqSTRd;
vqSTRdr=vqSTRd;

vqSTRdl(:,1:20,:)=0;
vqSTRdr(:,21:41,:)=0;

vqHPFl=vqHPF;
vqHPFr=vqHPF;

vqHPFl(:,1:20,:)=0;
vqHPFr(:,21:41,:)=0;

vqHYl=vqHY;
vqHYr=vqHY;

vqHYl(:,1:20,:)=0;
vqHYr(:,21:41,:)=0;

vqDORsml=vqDORsm;
vqDORsmr=vqDORsm;

vqDORsml(:,1:20,:)=0;
vqDORsmr(:,21:41,:)=0;

vqDORpml=vqDORpm;
vqDORpmr=vqDORpm;

vqDORpml(:,1:20,:)=0;
vqDORpmr(:,21:41,:)=0;

vqRSPl=vqRSP;
vqRSPr=vqRSP;

vqRSPl(:,1:20,:)=0;
vqRSPr(:,21:41,:)=0;

vqMDl=vqMD;
vqMDr=vqMD;

vqMDl(:,1:20,:)=0;
vqMDr(:,21:41,:)=0;


%%

%get boundary
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


%% 20230830 use atlas => change to white

%vqwtavatlas=vqwtav;
vqwtavatlas=vqwtav.*(~bdMD).*(~bdTRN).*(~bdACC).*(~bdPL).*(~bdIL).*(~bdRSP).*(~bdSTRd).*(~bdSTRv).*(~bdRE).*(~bdMO).*(~bdHPF).*(~bdHY).*(~bdDORsm).*(~bdDORpm);

%vqhetavatlas=vqhetav;
vqhetavatlas=vqhetav.*(~bdMD).*(~bdTRN).*(~bdACC).*(~bdPL).*(~bdIL).*(~bdRSP).*(~bdSTRd).*(~bdSTRv).*(~bdRE).*(~bdMO).*(~bdHPF).*(~bdHY).*(~bdDORsm).*(~bdDORpm);


%%
figure()

for i=27:35%19:1:43%21:30

    subplot(4,9,(i-27)/1+1)
    imagesc(permute(reshape(vqwtav(:,i,:),size(vqwtav,1),[]),[2,1]),[-0.2 0.5]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('WT')
    end
    
    subplot(4,9,(i-27)/1+1+9)
    imagesc(permute(reshape(vqhetav(:,i,:),size(vqhetav,1),[]),[2,1]),[-0.2 0.5]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('HET')
    end

    subplot(4,9,(i-27)/1+1+9*2)
    %imagesc(permute(reshape(vqwtavatlas(:,i,:),size(vqwtav,1),[]),[2,1]),[0.1 0.8]);
    
    imagesc(permute(reshape(vqwtav(:,i,:)+0.3*vqMD(:,i,:)+0.3*vqACC(:,i,:)+0.3*vqPL(:,i,:)+0.3*vqIL(:,i,:)+0.3*vqRSP(:,i,:)+0.3*vqSTRv(:,i,:)+0.3*vqSTRd(:,i,:)+0.3*vqRE(:,i,:),size(vqwtav,1),[]),[2,1]),[0.1 0.8]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('WT')
    end
    
    subplot(4,9,(i-27)/1+1+9*3)
    imagesc(permute(reshape(vqhetav(:,i,:)+0.3*vqMD(:,i,:)+0.3*vqACC(:,i,:)+0.3*vqPL(:,i,:)+0.3*vqIL(:,i,:)+0.3*vqRSP(:,i,:)+0.3*vqSTRv(:,i,:)+0.3*vqSTRd(:,i,:)+0.3*vqRE(:,i,:),size(vqhetav,1),[]),[2,1]),[0.1 0.8]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==11
    title('HET')
    end


end




%%
figure()

colormap("pink")

for i=1:3:41% 26:36%19:1:43%21:30

    subplot(4,14,(i-1)/3+1)
    %imagesc(permute(reshape(vqwtav(:,i,:),size(vqwtav,1),[]),[2,1]),[0.2 0.38]);
    imagesc(imadjust(permute(reshape(vqwtav(:,i,:),size(vqwtav,1),[]),[2,1]),[0 1],[0 1],0.7),[0.4 0.9]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==19
    title({'WT','Right'})
    elseif i==43
        title({'','Left'})
    else

        %title(strcat(num2str((i-31)*0.02),'mm'));
    end
    
    subplot(4,14,(i-1)/3+1+14)
    %imagesc(permute(reshape(vqhetav(:,i,:),size(vqhetav,1),[]),[2,1]),[0.2 0.38]);
    imagesc(imadjust(permute(reshape(vqhetav(:,i,:),size(vqhetav,1),[]),[2,1]),[0 1],[0 1],0.7),[0.4 0.9]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==19
    title('HET')
    end

    subplot(4,14,(i-1)/3+1+14*2)
    imagesc(imadjust(permute(reshape(vqwtavatlas(:,i,:),size(vqwtav,1),[]),[2,1]),[0 1],[0 1],0.7),[0.4 0.9]);
    %imagesc(permute(reshape(vqwtav(:,i,:)+1*bdMD(:,i,:)+1*bdACC(:,i,:)+1*bdPL(:,i,:)+1*bdIL(:,i,:)+1*bdRSP(:,i,:)+1*bdSTR(:,i,:),size(vqwtav,1),[]),[2,1]),[0.2 0.4]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==19
    title('WT')
    end
    
    subplot(4,14,(i-1)/3+1+14*3)
    imagesc(imadjust(permute(reshape(vqhetavatlas(:,i,:),size(vqwtav,1),[]),[2,1]),[0 1],[0 1],0.7),[0.4 0.9]);
    %imagesc(permute(reshape(vqhetav(:,i,:)+1*bdMD(:,i,:)+1*bdACC(:,i,:)+1*bdPL(:,i,:)+1*bdIL(:,i,:)+1*bdRSP(:,i,:)+1*bdSTR(:,i,:),size(vqhetav,1),[]),[2,1]),[0.2 0.4]);
    %xlim([10 53])
    %ylim([10 111])
    set(gca, 'YDir','normal');
    if i==19
    title('HET')
    end


end
%% subtraction
figure()

%colormap("pink")

for i=19:3:41% 26:36%19:1:43%21:30

    subplot(4,9,(i-19)/3+1)
    imagesc(permute(reshape(vqwtav(:,i,:)-vqhetav(:,i,:),size(vqwtav,1),[]),[2,1]),[0 0.18]);
    xlim([10 53])
    ylim([10 111])
    set(gca, 'YDir','normal');
    if i==19
    title({'WT-HET','Right'})
    elseif i==43
        title({'','Left'})
    else

        %title(strcat(num2str((i-31)*0.02),'mm'));
    end
    
    subplot(4,9,(i-19)/3+1+9)
    imagesc(permute(reshape(vqhetav(:,i,:)-vqwtav(:,i,:),size(vqhetav,1),[]),[2,1]),[0 0.18]);
    xlim([10 53])
    ylim([10 111])
    set(gca, 'YDir','normal');
    if i==19
    title('HET-WT')
    end

    subplot(4,9,(i-19)/3+1+9*2)
    %imagesc(permute(reshape(vqwtavatlas(:,i,:),size(vqwtav,1),[]),[2,1]),[0.2 0.38]);
    imagesc(permute(reshape(vqwtav(:,i,:)-vqhetav(:,i,:)+1*bdMD(:,i,:)+1*bdACC(:,i,:)+1*bdPL(:,i,:)+1*bdIL(:,i,:)+1*bdRSP(:,i,:)+1*bdSTR(:,i,:)+1*bdRE(:,i,:),size(vqwtav,1),[]),[2,1]),[0 0.2]);
    xlim([10 53])
    ylim([10 111])
    set(gca, 'YDir','normal');
    if i==19
    title('WT-HET')
    end
    
    subplot(4,9,(i-19)/3+1+9*3)
    imagesc(permute(reshape(vqhetavatlas(:,i,:),size(vqwtav,1),[]),[2,1]),[0.2 0.38]);
    imagesc(permute(reshape(vqhetav(:,i,:)-vqwtav(:,i,:)+1*bdMD(:,i,:)+1*bdACC(:,i,:)+1*bdPL(:,i,:)+1*bdIL(:,i,:)+1*bdRSP(:,i,:)+1*bdSTR(:,i,:)+1*bdRE(:,i,:),size(vqhetav,1),[]),[2,1]),[0 0.2]);
    xlim([10 53])
    ylim([10 111])
    set(gca, 'YDir','normal');
    if i==19
    title('HET-WT')
    end


end

%% AVERAGE SELECTED REGIONS

%define anteriorACC (where PL appears)
% vqACCa=vqACC;
% vqACCa(1:40,:,:)=0;

vqwtMOlav=zeros(length(wtid),1);
vqwtMOrav=zeros(length(wtid),1);
vqwtRSPlav=zeros(length(wtid),1);
vqwtRSPrav=zeros(length(wtid),1);
vqwtSTRdlav=zeros(length(wtid),1);
vqwtSTRdrav=zeros(length(wtid),1);
vqwtSTRvlav=zeros(length(wtid),1);
vqwtSTRvrav=zeros(length(wtid),1);
vqwtHPFlav=zeros(length(wtid),1);
vqwtHPFrav=zeros(length(wtid),1);
vqwtHYlav=zeros(length(wtid),1);
vqwtHYrav=zeros(length(wtid),1);
vqwtDORsmlav=zeros(length(wtid),1);
vqwtDORsmrav=zeros(length(wtid),1);
vqwtDORpmlav=zeros(length(wtid),1);
vqwtDORpmrav=zeros(length(wtid),1);

vqwtMDrav=zeros(length(wtid),1);
vqwtPLrav=zeros(length(wtid),1);
vqwtILrav=zeros(length(wtid),1);
vqwtACCrav=zeros(length(wtid),1);

vqwtTHrav=zeros(length(wtid),1);
vqwtPFCrav=zeros(length(wtid),1);
vqwtSTRrav=zeros(length(wtid),1);

vqhetMOlav=zeros(length(hetid),1);
vqhetMOrav=zeros(length(hetid),1);
vqhetRSPlav=zeros(length(hetid),1);
vqhetRSPrav=zeros(length(hetid),1);
vqhetSTRdlav=zeros(length(hetid),1);
vqhetSTRdrav=zeros(length(hetid),1);
vqhetSTRvlav=zeros(length(hetid),1);
vqhetSTRvrav=zeros(length(hetid),1);
vqhetHPFlav=zeros(length(hetid),1);
vqhetHPFrav=zeros(length(hetid),1);
vqhetHYlav=zeros(length(hetid),1);
vqhetHYrav=zeros(length(hetid),1);
vqhetDORsmlav=zeros(length(hetid),1);
vqhetDORsmrav=zeros(length(hetid),1);
vqhetDORpmlav=zeros(length(hetid),1);
vqhetDORpmrav=zeros(length(hetid),1);

vqhetMDrav=zeros(length(hetid),1);
vqhetPLrav=zeros(length(hetid),1);
vqhetILrav=zeros(length(hetid),1);
vqhetACCrav=zeros(length(hetid),1);

vqhetTHrav=zeros(length(hetid),1);
vqhetPFCrav=zeros(length(hetid),1);
vqhetSTRrav=zeros(length(hetid),1);

for wti=1:length(wtid_kcc)
    vqwtMOlav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqMOl)))/sum(sum(sum(vqMOl)));
    vqwtMOrav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqMOr)))/sum(sum(sum(vqMOr)));
    vqwtRSPlav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqRSPl)))/sum(sum(sum(vqRSPl)));
    vqwtRSPrav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqRSPr)))/sum(sum(sum(vqRSPr)));
    vqwtSTRvlav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqSTRvl)))/sum(sum(sum(vqSTRvl)));
    vqwtSTRvrav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqSTRvr)))/sum(sum(sum(vqSTRvr)));
    vqwtSTRdlav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqSTRdl)))/sum(sum(sum(vqSTRdl)));
    vqwtSTRdrav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqSTRdr)))/sum(sum(sum(vqSTRdr)));
    
    vqwtHPFlav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqHPFl)))/sum(sum(sum(vqHPFl)));
    vqwtHPFrav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqHPFr)))/sum(sum(sum(vqHPFr)));
    vqwtHYlav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqHYl)))/sum(sum(sum(vqHYl)));
    vqwtHYrav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqHYr)))/sum(sum(sum(vqHYr)));
    
    vqwtDORsmlav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqDORsml)))/sum(sum(sum(vqDORsml)));
    vqwtDORsmrav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqDORsmr)))/sum(sum(sum(vqDORsmr)));
    vqwtDORpmlav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqDORpml)))/sum(sum(sum(vqDORpml)));
    vqwtDORpmrav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqDORpmr)))/sum(sum(sum(vqDORpmr)));

    vqwtMDrav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqMDr)))/sum(sum(sum(vqMDr)));
    vqwtPLrav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqPLr)))/sum(sum(sum(vqPLr)));
    vqwtILrav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqILr)))/sum(sum(sum(vqILr)));
    vqwtACCrav(wti)=sum(sum(sum(vqwt(:,:,:,wti).*vqACCr)))/sum(sum(sum(vqACCr)));

    vqwtTHrav(wti)=(sum(sum(sum(vqwt(:,:,:,wti).*vqDORpmr)))+sum(sum(sum(vqwt(:,:,:,wti).*vqDORsmr))))/(sum(sum(sum(vqDORpmr)))+sum(sum(sum(vqDORsmr))));
    vqwtPFCrav(wti)=(sum(sum(sum(vqwt(:,:,:,wti).*vqPLr)))+sum(sum(sum(vqwt(:,:,:,wti).*vqILr)))+sum(sum(sum(vqwt(:,:,:,wti).*vqACCr))))/(sum(sum(sum(vqPLr)))+sum(sum(sum(vqILr)))+sum(sum(sum(vqACCr))));
    vqwtSTRrav(wti)=(sum(sum(sum(vqwt(:,:,:,wti).*vqSTRvr)))+sum(sum(sum(vqwt(:,:,:,wti).*vqSTRdr))))/(sum(sum(sum(vqSTRvr)))+sum(sum(sum(vqSTRdr))));
    
end

for wti=1:length(hetid_kcc)
    vqhetMOlav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqMOl)))/sum(sum(sum(vqMOl)));
    vqhetMOrav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqMOr)))/sum(sum(sum(vqMOr)));
    vqhetRSPlav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqRSPl)))/sum(sum(sum(vqRSPl)));
    vqhetRSPrav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqRSPr)))/sum(sum(sum(vqRSPr)));
    vqhetSTRvlav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqSTRvl)))/sum(sum(sum(vqSTRvl)));
    vqhetSTRvrav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqSTRvr)))/sum(sum(sum(vqSTRvr)));
    vqhetSTRdlav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqSTRdl)))/sum(sum(sum(vqSTRdl)));
    vqhetSTRdrav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqSTRdr)))/sum(sum(sum(vqSTRdr)));
    
    vqhetHPFlav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqHPFl)))/sum(sum(sum(vqHPFl)));
    vqhetHPFrav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqHPFr)))/sum(sum(sum(vqHPFr)));
    vqhetHYlav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqHYl)))/sum(sum(sum(vqHYl)));
    vqhetHYrav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqHYr)))/sum(sum(sum(vqHYr)));
    
    vqhetDORsmlav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqDORsml)))/sum(sum(sum(vqDORsml)));
    vqhetDORsmrav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqDORsmr)))/sum(sum(sum(vqDORsmr)));
    vqhetDORpmlav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqDORpml)))/sum(sum(sum(vqDORpml)));
    vqhetDORpmrav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqDORpmr)))/sum(sum(sum(vqDORpmr)));
    
    vqhetMDrav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqMDr)))/sum(sum(sum(vqMDr)));
    vqhetPLrav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqPLr)))/sum(sum(sum(vqPLr)));
    vqhetILrav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqILr)))/sum(sum(sum(vqILr)));
    vqhetACCrav(wti)=sum(sum(sum(vqhet(:,:,:,wti).*vqACCr)))/sum(sum(sum(vqACCr)));
    
    vqhetTHrav(wti)=(sum(sum(sum(vqhet(:,:,:,wti).*vqDORpmr)))+sum(sum(sum(vqhet(:,:,:,wti).*vqDORsmr))))/(sum(sum(sum(vqDORpmr)))+sum(sum(sum(vqDORsmr))));
    vqhetPFCrav(wti)=(sum(sum(sum(vqhet(:,:,:,wti).*vqPLr)))+sum(sum(sum(vqhet(:,:,:,wti).*vqILr)))+sum(sum(sum(vqhet(:,:,:,wti).*vqACCr))))/(sum(sum(sum(vqPLr)))+sum(sum(sum(vqILr)))+sum(sum(sum(vqACCr))));
    vqhetSTRrav(wti)=(sum(sum(sum(vqhet(:,:,:,wti).*vqSTRvr)))+sum(sum(sum(vqhet(:,:,:,wti).*vqSTRdr))))/(sum(sum(sum(vqSTRvr)))+sum(sum(sum(vqSTRdr))));
end

%%
vqwtav_all=[vqwtMOlav,vqwtMOrav,vqwtRSPlav,vqwtRSPrav,vqwtSTRvlav,vqwtSTRvrav,vqwtSTRdlav,vqwtSTRdrav,vqwtHPFlav,vqwtHPFrav,vqwtHYlav,vqwtHYrav,vqwtDORsmlav,vqwtDORsmrav,vqwtDORpmlav,vqwtDORpmrav,vqwtTHrav, vqwtPFCrav, vqwtSTRrav, vqwtMDrav, vqwtPLrav, vqwtILrav, vqwtACCrav];
vqhetav_all=[vqhetMOlav,vqhetMOrav,vqhetRSPlav,vqhetRSPrav,vqhetSTRvlav,vqhetSTRvrav,vqhetSTRdlav,vqhetSTRdrav,vqhetHPFlav,vqhetHPFrav,vqhetHYlav,vqhetHYrav,vqhetDORsmlav,vqhetDORsmrav,vqhetDORpmlav,vqhetDORpmrav,vqhetTHrav, vqhetPFCrav, vqhetSTRrav, vqhetMDrav, vqhetPLrav, vqhetILrav, vqhetACCrav];

%vqwtav_all=[vqwtMOrav,vqwtRSPrav,vqwtSTRvrav,vqwtSTRdrav,vqwtHPFrav,vqwtHYrav,vqwtDORsmrav,vqwtDORpmrav,vqwtTHrav, vqwtPFCrav, vqwtSTRrav, vqwtMDrav];
%vqhetav_all=[vqhetMOrav,vqhetRSPrav,vqhetSTRvrav,vqhetSTRdrav,vqhetHPFrav,vqhetHYrav,vqhetDORsmrav,vqhetDORpmrav,vqhetTHrav, vqhetPFCrav, vqhetSTRrav, vqhetMDrav];

vqav_all=zeros(numel(Skcc),size(vqwtav_all,2));
vqav_all([wtid_kcc],:)=vqwtav_all;
vqav_all([hetid_kcc],:)=vqhetav_all;


%%
titlelabel=({'MOr','RSPr','STRvr','STRdr','HPFr','HYr','DORsmr','DORpmr','THr','PFCr','STRr','MDr'});
titlelabel=({'','MOr','','RSPr','','STRvr','','STRdr','','HPFr','','HYr','','DORsmr','','DORpmr','THr','PFCr','STRr','MDr','PLr','ILr','ACCr'});

for i=1:length(titlelabel)
    subplot(1,length(titlelabel),i)
    %boxplot([vqwtav_all(:,i),vqhetav_all(:,i)]);hold on;
    scatter(ones(size(vqwtav_all,1),1),vqwtav_all(:,i),'green');hold on;
    scatter(2*ones(size(vqhetav_all,1),1),vqhetav_all(:,i),'red');
    
%     scatter(ones(2,1),vqwtav_all(1:2,i),'green','o','filled');hold on;
%     scatter(2*ones(3,1),vqhetav_all([1,2,9],i),'red','o','filled');
    
    errorbar([1 2],[mean(vqwtav_all(:,i)) mean(vqhetav_all(:,i))],[std(vqwtav_all(:,i))/sqrt(wtcount) std(vqhetav_all(:,i))/sqrt(hetcount)])
    
    

    title(titlelabel{i});
    ylabel('local connectivity');
    legend({'wt','het'},'Location','southoutside');
    grid on;
    xticks([1:2]);
    yticks([-0.5:0.1:1]);
    ylim([-1 1])
    xlim([0.8 2.2])

    [~,pval]=ttest2(vqwtav_all(:,i),vqhetav_all(:,i));
    text(1.3,0.4,strcat('p=',num2str(pval)));

end
