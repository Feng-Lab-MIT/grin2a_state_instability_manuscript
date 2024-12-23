
load('\\Fenglab03\GRIN2A_cognitive_project\10272024 behavior data_revision\Sum_Data_All.mat');
%load('\\Fenglab03\GRIN2A_cognitive_project\10272024 behavior data_revision\SSFO_revision.mat');
%20221206 add a function to flaten the sel cells 


%%

WTselIDs={'756','983','984','M0131','M0110','2863','M0109','M0123','M0122'};
Grin2aselIDs={'20326','TT01','TT08','TT10','20325','2041','74826'};%{'204110','20326','20416','TT01','TT08','TT10'};
MDselIDs={'VG1','984','983','03074','M0131','744','756'};%{'VG1','984','983','03074','M0131'};
PLselIDs={'03022','03074','M0131','M0110'};
SSFOselIDs={'94124','94126','74826','20326','94128'};%{'94124','94126','74826','20326'};

grin2a_sub=structure_data_into_selanimalID(grin2a,Grin2aselIDs);
WT_sub=structure_data_into_selanimalID(WT,WTselIDs);
PL_ON_sub=structure_data_into_selanimalID(PL_ON,PLselIDs);
MD_ON_sub=structure_data_into_selanimalID(MD_ON,MDselIDs);
SSFO_ON_sub=structure_data_into_selanimalID(SSFO_ON_1,SSFOselIDs);
PL_OFF_sub=structure_data_into_selanimalID(PL_OFF,PLselIDs);
MD_OFF_sub=structure_data_into_selanimalID(MD_OFF,MDselIDs);
SSFO_OFF_sub=structure_data_into_selanimalID(SSFO_OFF_1,SSFOselIDs);


%%
save('behaviordata20241027.mat');

%%
%20221206 add a function to flaten the sel cells 
WT_sub_allcells=getallselcell(WT_sub,WTselIDs);
grin2a_sub_allcells=getallselcell(grin2a_sub,Grin2aselIDs);
SSFO_ONsub_allcells=getallselcell(SSFO_ON_sub,SSFOselIDs);
SSFO_OFFsub_allcells=getallselcell(SSFO_OFF_sub,SSFOselIDs);


%% need to plot # of switchs
%%
grin2a_sub=structure_data_into_animalID(grin2a);
WT_sub=structure_data_into_animalID(WT);
PL_ON_sub=structure_data_into_animalID(PL_ON);
MD_ON_sub=structure_data_into_animalID(MD_ON);
%SSFO_ON_sub=structure_data_into_animalID(SSFO_ON);
PL_OFF_sub=structure_data_into_animalID(PL_OFF);
MD_OFF_sub=structure_data_into_animalID(MD_OFF);
%SSFO_OFF_sub=structure_data_into_animalID(SSFO_OFF);

SSFO_ON_1_sub=structure_data_into_animalID(SSFO_ON_1);
SSFO_OFF_1_sub=structure_data_into_animalID(SSFO_OFF_1);

%%
% f=fieldnames(SSFO_OFF_1_sub);
% for i=1:length(f)
%     SSFO_OFF_sub.(f{i})=SSFO_OFF_1_sub.(f{i})
%     
% end
% 
% f=fieldnames(SSFO_ON_1_sub);
% for i=1:length(f)
%     SSFO_ON_sub.(f{i})=SSFO_ON_1_sub.(f{i})
%     
% end
%%
clearvars f i



%%

function datasub=structure_data_into_animalID(datacell)

animalIDWT=cellfun(@(x)x.ID,datacell,'UniformOutput',false);

unianimalIDWT=unique(animalIDWT);

for i=1:length(unianimalIDWT)
    datasub.(strcat('n',unianimalIDWT{i}))={datacell{strcmp(animalIDWT,unianimalIDWT{i})}};
end

end

function datasub=structure_data_into_selanimalID(datacell,selID)

animalIDWT=cellfun(@(x)x.ID,datacell,'UniformOutput',false);

%unianimalIDWT=unique(animalIDWT);

for i=1:length(selID)
    datasub.(strcat('n',selID{i}))={datacell{strcmp(animalIDWT,selID{i})}};
end

end


function An_sub_allcells=getallselcell(An_sub,selIDs)

for i=1:length(selIDs)
    if i==1
        An_sub_allcells=An_sub.(strcat('n',selIDs{i}));
    else
        An_sub_allcells={An_sub_allcells{:},An_sub.(strcat('n',selIDs{i})){:}}
    end
end

end
%%
% animalID=cellfun(@(x)x.ID,grin2a,'UniformOutput',false);
% 
% unianimalID=unique(animalID);
% 
% for i=1:length(unianimalID)
%     grin2a_sub.(strcat('n',unianimalID{i}))={grin2a{strcmp(animalID,unianimalID{i})}};
% end