%% summary plot with bootstramp so that each animal has equal contribution to the data
%%
% prepare data 
close all
clear all

load('Sum_Data_All.mat');
%grin2a=[grin2a,SSFO_grin2a_1];
n_run=100;

%% run WT v.s. Grin2a 
figure;title('WT v.s. Grin2a')
data=WT;
 IDs={'756','983','984','M0131','M0110','2863','M0109','M0123','M0122','M0119','M0121','VG1'};
%IDs={'756','983','984','M0131','M0110','2863'};
Y_WT_all=[];
M_WT_all=[];
M_BL_WT_all=[];
optimality_WT_all=[];
pswitch_WT_all=[];
ID_WT=[];
for j=1:length(IDs)
    %get current data input: IDs and data output: Current
Current=getCurrent(IDs{j},data);  

% get block length
BL=getBlockLength(Current);
M_BL_WT=mean(BL,'omitnan');
optimality=getOptimality(Current);
M_optimality_WT=mean(optimality,'omitnan');
% get policy
[Y_WT,M_WT]=GetMatrix(Current);


X=1:1:54;
Y_WT=Y_WT(:,1:54);
% bootstramp

% Y_WT=bootstrap(X,Y_WT,n_run);
Y_WT=smoothdata(Y_WT,2,'gaussian',2);
% plot
S_WT=estimateSEM(Y_WT);
M_WT=median(Y_WT,'omitnan');
subplot(4,9,j)
error_area(X,M_WT,S_WT,[0.5,0.5,0.5],0.5,7)
xlim([0 54])
ylim([0 1])

for i=1:size(Y_WT,1)
    sampleStart_WT(i)=min(find(Y_WT(i,:)<1));
    a=find(Y_WT(i,:)==0);
    if isempty(a)
        sampleEnd_WT(i)=nan;
    else
    sampleEnd_WT(i)=min(find(Y_WT(i,:)==0));
    end
end
M_sampleStart_WT(j)=mean(sampleStart_WT);
M_sampleEnd_WT(j)=mean(sampleEnd_WT);

Y_WT_all=[Y_WT_all;Y_WT];

M_WT_all=[M_WT_all; M_WT];

M_BL_WT_all(j)=M_BL_WT;
M_optimality_WT_all(j)=M_optimality_WT;


end



%%
data=grin2a;
IDs={'20326','TT01','TT08','TT10','20325','2041','74826'};% 2041 incluce 204110 and 20416
M_grin2a_all=[];
Y_grin2a_all=[];
% pswitch_grin2a_all=[];
M_BL_grin2a_all=[];
optimality_grin2a_all=[];

ID_grin2a=[];
for j=1:length(IDs)
Current=getCurrent(IDs{j},data);

% get block length
BL=getBlockLength(Current);
M_BL_grin2a=mean(BL,'omitnan');

optimality=getOptimality(Current);
M_optimality_grin2a=mean(optimality,'omitnan');
a=IDs{j}
subplot(4,9,j+18)
[Y_grin2a,M_grin2a]=GetMatrix(Current);
X=1:1:54;
Y_grin2a=Y_grin2a(:,1:54);
%BOOTSTRAP

% Y_grin2a=bootstrap(X,Y_grin2a,n_run);

Y_grin2a=smoothdata(Y_grin2a,2,'gaussian',2);
S_grin2a=estimateSEM(Y_grin2a);
M_grin2a=median(Y_grin2a,'omitnan');
error_area(X,M_grin2a,S_grin2a,[1,0.3,0.3],0.5,7)
xlim([0 54])
ylim([0 1])

for i=1:size(Y_grin2a,1)
    sampleStart_grin2a(i)=min(find(Y_grin2a(i,:)<1));
    a=find(Y_grin2a(i,:)==0);
    if isempty(a)
        sampleEnd_grin2a(i)=nan;
    else
    sampleEnd_grin2a(i)=min(find(Y_grin2a(i,:)==0));
    end
end
M_sampleStart_grin2a(j)=mean(sampleStart_grin2a);
M_sampleEnd_grin2a(j)=mean(sampleEnd_grin2a,'omitnan');

Y_grin2a_all=[Y_grin2a_all;Y_grin2a];

M_grin2a_all=[M_grin2a_all; M_grin2a];
% % get switch 
% [Y_pswitch]=getswitch(Current);
% a=size(Y_pswitch);
% ID_grin2a(j)=a(1);
% subplot(4,9,j+27);
% X=1:1:54;
% Y_pswitch=Y_pswitch(:,1:54);
% Y_pswitch=smoothdata(Y_pswitch,2,'gaussian',4);
% S_grin2a=estimateSEM(Y_pswitch);
% M_grin2a=mean(Y_pswitch,'omitnan');
% error_area(X,M_grin2a,S_grin2a,[1,0.3,0.3],0.5,7)
% xlim([0 54])
% ylim([0 0.7])
% 
% 
% % Y_pswitch=bootstrap(X,Y_pswitch,n_run);
% pswitch_grin2a_all=[pswitch_grin2a_all;Y_pswitch];

M_BL_grin2a_all(j)=M_BL_grin2a;
M_optimality_grin2a_all(j)=M_optimality_grin2a;
end


figure; 
M_M_WT_all=median(M_WT_all,'omitnan');
S_M_WT=estimateSEM(M_WT_all);
error_area(X,M_M_WT_all, S_M_WT, [0.5,0.5,0.5], 0.5,7)
xlim([0 54])
ylim([0 1])

hold on;

M_M_grin2a_all=median(M_grin2a_all,'omitnan');
S_M_grin2a=estimateSEM(M_grin2a_all);
error_area(X,M_M_grin2a_all, S_M_grin2a, [1,0.3,0.3], 0.5,7)

% %%
% figure; title('MD')
% 
% data=MD_OFF;
% IDs={'VG1','984','983','03074','M0131','744','756'};
% Y_MD_OFF_all=[];
% pswitch_MD_OFF_all=[];
% M_BL_MD_OFF_all=[];
% optimality_MD_OFF_all=[];
% 
% for j=1:length(IDs)
%     %get current data input: IDs and data output: Current
% Current=getCurrent(IDs{j},data);   
% % get block length
% BL=getBlockLength(Current);
% M_BL_MD_OFF=mean(BL,'omitnan');
% optimality=getOptimality(Current);
% M_optimality_MD_OFF=mean(optimality,'omitnan');
% 
% subplot(4,9,j)
% [Y_MD_OFF,M_MD_OFF]=GetMatrix(Current);
% X=1:1:54;
% Y_MD_OFF=Y_MD_OFF(:,1:54);
% % bootstramp
% 
% % Y_MD_OFF=bootstrap(X,Y_MD_OFF,n_run);
% Y_MD_OFF=smoothdata(Y_MD_OFF,2,'gaussian',2);
% % plot
% S_MD_OFF=estimateSEM(Y_MD_OFF);
% M_MD_OFF=median(Y_MD_OFF,'omitnan');
% error_area(X,M_MD_OFF,S_MD_OFF,[0.5,0.5,0.5],0.5,7)
% xlim([0 54])
% ylim([0 1])
% Y_MD_OFF_all=[Y_MD_OFF_all;Y_MD_OFF];
% M_MD_OFF_all=[M_MD_OFF_all; M_MD_OFF];
% 
% % % get switch 
% % [Y_pswitch]=getswitch(Current);
% % 
% % subplot(4,9,j+9);
% % X=1:1:54;
% % Y_pswitch=Y_pswitch(:,1:54);
% % Y_pswitch=smoothdata(Y_pswitch,2,'gaussian',4);
% % S_MD_OFF=estimateSEM(Y_pswitch);
% % M_MD_OFF=median(Y_pswitch,'omitnan');
% % error_area(X,M_MD_OFF,S_MD_OFF,[0.5,0.5,0.5],0.5,7)
% % xlim([0 54])
% % ylim([0 0.7])
% % 
% % 
% % %Y_pswitch=bootstrap(X,Y_pswitch,n_run);
% % pswitch_MD_OFF_all=[pswitch_MD_OFF_all;Y_pswitch];
% 
% M_BL_MD_OFF_all(j)=M_BL_MD_OFF;
% M_optimality_MD_OFF_all(j)=M_optimality_MD_OFF;
% 
% end
% 
% data=MD_ON;
% 
% Y_MD_ON_all=[];
% pswitch_MD_ON_all=[];
% M_BL_MD_ON_all=[];
% optimality_MD_ON_all=[];
% 
% for j=1:length(IDs)
% Current=getCurrent(IDs{j},data);   
% 
% % get block length
% BL=getBlockLength(Current);
% M_BL_MD_ON=mean(BL,'omitnan');
% optimality=getOptimality(Current);
% M_optimality_MD_ON=mean(optimality,'omitnan');
% 
% subplot(4,9,j+18)
% [Y_MD_ON,M_MD_ON]=GetMatrix(Current);
% X=1:1:54;
% Y_MD_ON=Y_MD_ON(:,1:54);
% % BOOTSTRAP
% % Y_MD_ON=bootstrap(X,Y_MD_ON,n_run);
% %
% Y_MD_ON=smoothdata(Y_MD_ON,2,'gaussian',2);
% S_MD_ON=estimateSEM(Y_MD_ON);
% M_MD_ON=median(Y_MD_ON,'omitnan');
% plot(X,M_MD_ON)
% %error_area(X,M_MD_ON,S_MD_ON,[1,0.3,0.3],0.5,7)
% xlim([0 54])
% ylim([0 1])
% 
% Y_MD_ON_all=[Y_MD_ON_all;Y_MD_ON];
% M_MD_ON_all=[M_MD_ON_all; M_MD_ON];
% % % get switch 
% % [Y_pswitch]=getswitch(Current);
% % 
% % subplot(4,9,j+27);
% % X=1:1:54;
% % Y_pswitch=Y_pswitch(:,1:54);
% % Y_pswitch=smoothdata(Y_pswitch,2,'gaussian',4);
% % S_MD_ON=estimateSEM(Y_pswitch);
% % 
% % M_MD_ON=median(Y_pswitch,'omitnan');
% % %error_area(X,M_MD_ON,S_MD_ON,[1,0.3,0.3],0.5,7)
% % xlim([0 54])
% % ylim([0 0.7])
% % 
% % 
% % %Y_pswitch=bootstrap(X,Y_pswitch,n_run);
% % pswitch_MD_ON_all=[pswitch_MD_ON_all;Y_pswitch];
% 
% M_BL_MD_ON_all(j)=M_BL_MD_ON;
% M_optimality_MD_ON_all(j)=M_optimality_MD_ON;
% 
% end
% % 
% 
% 
% %%
% figure; title('PL')
% data=PL_OFF;
% IDs={'03022','03074','M0131','M0110'};
% Y_PL_OFF_all=[];
% pswitch_PL_OFF_all=[];
% for j=1:length(IDs)
%     %%get current data input: IDs and data output: Current
% Current=getCurrent(IDs{j},data);   
% subplot(4,9,j)
% [Y_PL_OFF,M_PL_OFF]=GetMatrix(Current); 
% X=1:1:54;
% Y_PL_OFF=Y_PL_OFF(:,1:54);
% %bootstramp
% 
% % Y_PL_OFF=bootstrap(X,Y_PL_OFF,n_run);
% Y_PL_OFF=smoothdata(Y_PL_OFF,2,'gaussian',2);
% %plot
% S_PL_OFF=estimateSEM(Y_PL_OFF);
% M_PL_OFF=median(Y_PL_OFF,'omitnan');
% error_area(X,M_PL_OFF,S_PL_OFF,[0.5,0.5,0.5],0.5,7)
% xlim([0 54])
% ylim([0 1])
% Y_PL_OFF_all=[Y_PL_OFF_all;Y_PL_OFF];
% 
% % % get switch 
% % [Y_pswitch]=getswitch(Current);
% % 
% % subplot(4,9,j+9);
% % X=1:1:54;
% % Y_pswitch=Y_pswitch(:,1:54);
% % Y_pswitch=smoothdata(Y_pswitch,2,'gaussian',4);
% % S_PL_OFF=estimateSEM(Y_pswitch);
% % M_PL_OFF=median(Y_pswitch,'omitnan');
% % error_area(X,M_PL_OFF,S_PL_OFF,[0.5,0.5,0.5],0.5,7)
% % xlim([0 54])
% % ylim([0 0.7])
% % 
% % 
% % %Y_pswitch=bootstrap(X,Y_pswitch,n_run);
% % pswitch_PL_OFF_all=[pswitch_PL_OFF_all;Y_pswitch];
% 
% end
% 
% data=PL_ON;
% 
% Y_PL_ON_all=[];
% pswitch_PL_ON_all=[];
% for j=1:length(IDs)
% Current=getCurrent(IDs{j},data);   
% subplot(4,9,j+18)
% [Y_PL_ON,M_PL_ON]=GetMatrix(Current);
% X=1:1:54;
% Y_PL_ON=Y_PL_ON(:,1:54);
% % BOOTSTRAP
% 
% % Y_PL_ON=bootstrap(X,Y_PL_ON,n_run);
% %
% Y_PL_ON=smoothdata(Y_PL_ON,2,'gaussian',2);
% S_PL_ON=estimateSEM(Y_PL_ON);
% M_PL_ON=median(Y_PL_ON,'omitnan');
% error_area(X,M_PL_ON,S_PL_ON,[1,0.3,0.3],0.5,7)
% xlim([0 54])
% ylim([0 1])
% 
% Y_PL_ON_all=[Y_PL_ON_all;Y_PL_ON];
% 
% % % get switch 
% % [Y_pswitch]=getswitch(Current);
% % 
% % subplot(4,9,j+27);
% % X=1:1:54;
% % Y_pswitch=Y_pswitch(:,1:54);
% % Y_pswitch=smoothdata(Y_pswitch,2,'gaussian',4);
% % S_PL_ON=estimateSEM(Y_pswitch);
% % M_PL_ON=median(Y_pswitch,'omitnan');
% % error_area(X,M_PL_ON,S_PL_ON,[1,0.3,0.3],0.5,7)
% % xlim([0 54])
% % ylim([0 0.7])
% % 
% % 
% % %Y_pswitch=bootstrap(X,Y_pswitch,n_run);
% % pswitch_PL_ON_all=[pswitch_PL_ON_all;Y_pswitch];
% end
% %%
% figure; title('SSFO')
% data=SSFO_OFF_1;
% IDs={'94124','94126','74826','20326'};
% Y_SSFO_OFF_1_all=[];
% pswitch_SSFO_OFF_all=[];
% for j=1:length(IDs)
%     %%get current data input: IDs and data output: Current
% Current=getCurrent(IDs{j},data);   
% subplot(4,9,j)
% [Y_SSFO_OFF_1,M_SSFO_OFF_1]=GetMatrix(Current);
% X=1:1:54;
% Y_SSFO_OFF_1=Y_SSFO_OFF_1(:,1:54);
% % bootstramp
% % Y_SSFO_OFF_1=bootstrap(X,Y_SSFO_OFF_1,n_run);
% Y_SSFO_OFF_1=smoothdata(Y_SSFO_OFF_1,2,'gaussian',2);
% % plot
% S_SSFO_OFF_1=estimateSEM(Y_SSFO_OFF_1);
% M_SSFO_OFF_1=median(Y_SSFO_OFF_1,'omitnan');
% error_area(X,M_SSFO_OFF_1,S_SSFO_OFF_1,[0.5,0.5,0.5],0.5,7)
% xlim([0 54])
% ylim([0 1])
% Y_SSFO_OFF_1_all=[Y_SSFO_OFF_1_all;Y_SSFO_OFF_1];
% % % get switch 
% % [Y_pswitch]=getswitch(Current);
% % 
% % subplot(4,9,j+9);
% % X=1:1:54;
% % Y_pswitch=Y_pswitch(:,1:54);
% % Y_pswitch=smoothdata(Y_pswitch,2,'gaussian',4);
% % S_SSFO_OFF=estimateSEM(Y_pswitch);
% % M_SSFO_OFF=median(Y_pswitch,'omitnan');
% % error_area(X,M_SSFO_OFF,S_SSFO_OFF,[0.5,0.5,0.5],0.5,7)
% % xlim([0 54])
% % ylim([0 0.7])
% % 
% % 
% % %Y_pswitch=bootstrap(X,Y_pswitch,n_run);
% % pswitch_SSFO_OFF_all=[pswitch_SSFO_OFF_all;Y_pswitch];
% 
% end
% 
% data=SSFO_ON_1;
% 
% Y_SSFO_ON_1_all=[];
% pswitch_SSFO_ON_all=[];
% for j=1:length(IDs)
% Current=getCurrent(IDs{j},data);   
% subplot(4,9,j+18)
% [Y_SSFO_ON_1,M_SSFO_ON_1]=GetMatrix(Current);
% X=1:1:54;
% Y_SSFO_ON_1=Y_SSFO_ON_1(:,1:54);
% % BOOTSTRAP
% % Y_SSFO_ON_1=bootstrap(X,Y_SSFO_ON_1,n_run);
% %
% Y_SSFO_ON_1=smoothdata(Y_SSFO_ON_1,2,'gaussian',2);
% S_SSFO_ON_1=estimateSEM(Y_SSFO_ON_1);
% M_SSFO_ON_1=median(Y_SSFO_ON_1,'omitnan');
% error_area(X,M_SSFO_ON_1,S_SSFO_ON_1,[1,0.3,0.3],0.5,7)
% xlim([0 54])
% ylim([0 1])
% 
% Y_SSFO_ON_1_all=[Y_SSFO_ON_1_all;Y_SSFO_ON_1];
% % % get switch 
% % [Y_pswitch]=getswitch(Current);
% % 
% % subplot(4,9,j+27);
% % X=1:1:54;
% % Y_pswitch=Y_pswitch(:,1:54);
% % Y_pswitch=smoothdata(Y_pswitch,2,'gaussian',4);
% % S_SSFO_ON=estimateSEM(Y_pswitch);
% % M_SSFO_ON=median(Y_pswitch,'omitnan');
% % error_area(X,M_SSFO_ON,S_SSFO_ON,[1,0.3,0.3],0.5,7)
% % xlim([0 54])
% % ylim([0 0.7])
% % 
% % 
% % %Y_pswitch=bootstrap(X,Y_pswitch,n_run);
% % pswitch_SSFO_ON_all=[pswitch_SSFO_ON_all;Y_pswitch];
% 
% end
% %
% %%
% figure;
% subplot(2,4,1)
% 
% S_WT_all=estimateSEM(Y_WT_all);
% M_WT_all=median(Y_WT_all,'omitnan');
% error_area(X,M_WT_all,S_WT_all,[0.3,0.3,0.3],0.5,7)
% xlim([0 54])
% ylim([0 1])
% hold on;
% S_grin2a_all=estimateSEM(Y_grin2a_all);
% M_grin2a_all=median(Y_grin2a_all,'omitnan');
% error_area(X,M_grin2a_all,S_grin2a_all,[1,0.3,0.3],0.5,7)
% title ('WT Grin2a')
% 
% % subplot(2,4,5)
% % 
% % S_pswitch_WT_all=estimateSEM(pswitch_WT_all);
% % M_pswitch_WT_all=mean(pswitch_WT_all,'omitnan');
% % error_area(X,M_pswitch_WT_all,S_pswitch_WT_all,[0.3,0.3,0.3],0.5,7)
% % xlim([0 54])
% % ylim([0 0.7])
% % hold on;
% % S_pswitch_grin2a_all=estimateSEM(pswitch_grin2a_all);
% % M_pswitch_grin2a_all=mean(pswitch_grin2a_all,'omitnan');
% % error_area(X,M_pswitch_grin2a_all,S_pswitch_grin2a_all,[1,0.3,0.3],0.5,7)
% % title ('WT Grin2a p(switch)')
% % 
% 
% %
% subplot(2,4,4)
% 
% S_SSFO_OFF_1_all=estimateSEM(Y_SSFO_OFF_1_all);
% M_SSFO_OFF_1_all=median(Y_SSFO_OFF_1_all,'omitnan');
% error_area(X,M_SSFO_OFF_1_all,S_SSFO_OFF_1_all,[0.3,0.3,0.3],0.5,7)
% xlim([0 54])
% ylim([0 1])
% hold on;
% S_SSFO_ON_1_all=estimateSEM(Y_SSFO_ON_1_all);
% M_SSFO_ON_1_all=median(Y_SSFO_ON_1_all,'omitnan');
% error_area(X,M_SSFO_ON_1_all,S_SSFO_ON_1_all,[1,0.3,0.3],0.5,7)
% title('SSFO_grin2a')
% 
% % subplot(2,4,8)
% % 
% % S_pswitch_SSFO_OFF_all=estimateSEM(pswitch_SSFO_OFF_all);
% % M_pswitch_SSFO_OFF_all=mean(pswitch_SSFO_OFF_all,'omitnan');
% % error_area(X,M_pswitch_SSFO_OFF_all,S_pswitch_SSFO_OFF_all,[0.3,0.3,0.3],0.5,7)
% % xlim([0 54])
% % ylim([0 0.7])
% % hold on;
% % S_pswitch_SSFO_ON_all=estimateSEM(pswitch_SSFO_ON_all);
% % M_pswitch_SSFO_ON_all=mean(pswitch_SSFO_ON_all,'omitnan');
% % error_area(X,M_pswitch_SSFO_ON_all,S_pswitch_SSFO_ON_all,[1,0.3,0.3],0.5,7)
% % 
% % title ('SSFO Grin2a p(switch)')
% 
% %%
% 
% subplot(2,4,3)
% 
% S_MD_OFF_all=estimateSEM(Y_MD_OFF_all);
% M_MD_OFF_all=median(Y_MD_OFF_all,'omitnan');
% error_area(X,M_MD_OFF_all,S_MD_OFF_all,[0.3,0.3,0.3],0.5,7)
% 
% hold on;
% S_MD_ON_all=estimateSEM(Y_MD_ON_all);
% M_MD_ON_all=median(Y_MD_ON_all,'omitnan');
% error_area(X,M_MD_ON_all,S_MD_ON_all,[1,0.3,0.3],0.5,7)
% title ('MD')
% 
% % subplot(2,4,7)
% % 
% % S_pswitch_MD_OFF_all=estimateSEM(pswitch_MD_OFF_all);
% % M_pswitch_MD_OFF_all=mean(pswitch_MD_OFF_all,'omitnan');
% % error_area(X,M_pswitch_MD_OFF_all,S_pswitch_MD_OFF_all,[0.3,0.3,0.3],0.5,7)
% % xlim([0 54])
% % ylim([0 0.7])
% % hold on;
% % S_pswitch_MD_ON_all=estimateSEM(pswitch_MD_ON_all);
% % M_pswitch_MD_ON_all=mean(pswitch_MD_ON_all,'omitnan');
% % error_area(X,M_pswitch_MD_ON_all,S_pswitch_MD_ON_all,[1,0.3,0.3],0.5,7)
% % 
% % 
% % 
% 
% subplot(2,4,2)
% 
% 
% S_PL_OFF_all=estimateSEM(Y_PL_OFF_all);
% M_PL_OFF_all=median(Y_PL_OFF_all,'omitnan');
% error_area(X,M_PL_OFF_all,S_PL_OFF_all,[0.3,0.3,0.3],0.5,7)
% 
% hold on;
% S_PL_ON_all=estimateSEM(Y_PL_ON_all);
% M_PL_ON_all=median(Y_PL_ON_all,'omitnan');
% error_area(X,M_PL_ON_all,S_PL_ON_all,[1,0.3,0.3],0.5,7)
% 
% title('PL')
% % 
% % subplot(2,4,6)
% % 
% % S_pswitch_PL_OFF_all=estimateSEM(pswitch_PL_OFF_all);
% % M_pswitch_PL_OFF_all=mean(pswitch_PL_OFF_all,'omitnan');
% % error_area(X,M_pswitch_PL_OFF_all,S_pswitch_PL_OFF_all,[0.3,0.3,0.3],0.5,7)
% % xlim([0 54])
% % ylim([0 0.7])
% % hold on;
% % S_pswitch_PL_ON_all=estimateSEM(pswitch_PL_ON_all);
% % M_pswitch_PL_ON_all=mean(pswitch_PL_ON_all,'omitnan');
% % error_area(X,M_pswitch_PL_ON_all,S_pswitch_PL_ON_all,[1,0.3,0.3],0.5,7)
% % 
% 
% % %% summary plot without bootstramp
% % % summary plot
% % 
% %  figure;
% %  subplot(1,4,1)
% %  plotY(MD_WT,'ON')
% %  hold on;
% %  plotY(MD_grin2a,'OFF')
% %  ylim([0 1])
% %  xlim([0 40])
% %  title('MD')
% % 
% %  subplot(1,4,2)
% %   plotY(PL_WT,'ON')
% %  hold on;
% %  plotY(PL_grin2a,'OFF')
% %  ylim([0 1])
% %  xlim([0 40])
% %  title('PFC')
% % 
% %  subplot(1,4,3)
% %  plotY(WT,'OFF')
% %  hold on;
% %  plotY(grin2a,'grin2a')
% %  ylim([0 1])
% %  xlim([0 40])
% %  title('WT v.s.grin2a')
% % 
% %   subplot(1,4,4)
% %  plotY(SSFO_grin2a,'OFF')
% %  hold on;
% %  plotY(SSFO_WT,'ON')
% %  ylim([0 1])
% %  xlim([0 40])
% %  title('grin2aSSFO ON v.s.OFF')
% % 
% % 
% % 
% %  figure;
% %  subplot(1,4,1)
% %  plotMean(MD_WT,'ON')
%  hold on;
%  plotMean(MD_grin2a,'OFF')
%  ylim([0 1])
%  title('MD')
% 
%  subplot(1,4,2)
%   plotMean(PL_WT,'ON')
%  hold on;
%  plotMean(PL_grin2a,'OFF')
%  ylim([0 1])
%  title('PFC')
% 
%  subplot(1,4,3)
%  plotMean(WT,'OFF')
%  hold on;
%  plotMean(grin2a,'grin2a')
%  ylim([0 1])
%  title('WT v.s.grin2a')
% 
%   subplot(1,4,4)
%  plotMean(SSFO_grin2a_1,'OFF')
%  hold on;
%  plotMean(SSFO_WT_1,'ON')
%  ylim([0 1])
%  title('grin2aSSFO ON v.s.OFF')



 %% functions
 
 
 
%plot MD_WT and PL_WT


function plotY(data,type)
[Y,M]=GetMatrix(data);

Y=smoothdata(Y,2,'gaussian',2);
    X=1:1:54;
    Y=Y(:,1:54);

S=estimateSEM(Y);
M=median(Y,'omitnan');

if contains(type,'ON')
 error_area(X,M,S,[0.7,0.7,0.1],0.5,7)
else
if contains(type,'grin2a')
    error_area(X,M,S,[0.7,0.1,0.1],0.5,7)
else
    error_area(X,M,S,[0.3,0.3,0.3],0.5,7)
end
end
end


function plotMean(data,type)
[Y,M]=GetMatrix(data);

Y=smoothdata(Y,2,'gaussian',2);
    X=1:1:54;
    Y=Y(:,1:54);

S=estimateSEM(Y);
M=mean(Y,'omitnan');

if contains(type,'ON')
 error_area(X,M,S,[0.7,0.7,0.1],0.5,7)
else
if contains(type,'grin2a')
    error_area(X,M,S,[0.7,0.1,0.1],0.5,7)
else
    error_area(X,M,S,[0.3,0.3,0.3],0.5,7)
end
end
end





function [Y_pHR_matrix,Y_model_matrix] = GetMatrix(block)



Y_pHR_matrix=[];
Y_model_matrix=[];
for i=1:length(block)
    if ~isempty(block)
    currentHRrequest=block{i}.HRrequest;
    currentHRchoice=block{i}.HRchoice;
    
[X_HRrequest,Y_pHR] = estimatePolicy(currentHRrequest,currentHRchoice);
[k,model,Y_model] = fitSigmoid(Y_pHR,X_HRrequest);

n=length(Y_pHR);

if testBlock(block{i})==1
    Y_pHR(n:54)=0;
else
    Y_pHR(n:54)=nan;
end

if length(Y_pHR)>54
    Y_pHR(55:end)=[];
end
% if Y_pHR(1)>=1
if isnan(Y_pHR(1))
    Y_pHR(1)=1;
end
%  Onset(i)=min(find(Y_pHR<1));  
    else
        Y_pHR=NaN(1,54);
        Y_model=NAN(1,54);
    end
 Y_pHR_matrix(i,:)=Y_pHR;
 Y_model_matrix(i,:)=Y_model;
 
 
end
end

 function [k,model,Y_model] = fitSigmoid(Y,X)
    % Input only Y: Y=HRchoice
    % Inputs inlcude Y and X: Y=Y_pHR; X=X_HRrequest
    warning off
    type = 'Normal';
    model = @(k,x)  k(3).*(1-cdf(type,x,k(1),k(2)))+(1-k(3));
    try
        if nargin==1
            if size(Y,2)>size(Y,1)
                Y = Y';
            else
            end
            k = lsqcurvefit(model,[length(Y)/2 3 1],1:length(Y),Y',[15 1 1],[50 5 1],optimset('Display','off'));
            Y_model = model(k,[1:100]');
        else
            if size(Y,2)>size(Y,1)
                Y = Y';
            else
            end
            if size(X,2)>size(X,1)
                X = X';
            else
            end
            k = lsqcurvefit(model,[18 3 1],X',Y',[10 1 1],[30 5 1],optimset('Display','off'));
            %3 parameters of initial [midpoint slope LRcommitment] while [10 0.5 0.5] and [30 5 1] are lower and upper boundaries, respectively.
            Y_model = model(k,[1:100]');
        end
    catch
        k = nan(1,3);
        Y_model = nan(1,100);
    end
 end
 
 function Y=bootstrap(X,Y,n_run)

x_HR_req_list=X;
n_animal = 1;
HR_LR_animal_bs_list = {}; % 1=HR, 0=LR.
for i_animal = 1
    HR_LR_bs_list_temp = NaN(n_run,length(x_HR_req_list));
    for i_HR_req = 1:length(x_HR_req_list)
        x_HR_req = x_HR_req_list(i_HR_req);
        HR_LR_bs_list_temp(:,i_HR_req) = randsample(Y(:,i_HR_req),n_run,true);
    end
end


Y=HR_LR_bs_list_temp;
 end
 
 function Current=getCurrent(ID,data)
 AnimalID=ID;
    eval(['WT' AnimalID '={}']);
    
    for i=1:length(data)
        if ~isempty(data{i})
    if contains(data{i}.ID,AnimalID)
      eval(['WT' AnimalID '{end+1}' '=data{i}']); % it worked!!!!
    end
        end
    end   
    eval(['Current=WT' AnimalID]);  
 end

 function [Y_pswitch]=getswitch(block);

Y_pswitch=zeros(length(block),200);

for i=1:length(block)

    if ~isempty(block)
    HRrequest=block{i}.HRrequest;
    LRchoice=block{i}.LRchoice;

    X_HRrequest = unique(HRrequest);
    Y_switch=diff(LRchoice);
    Y_switch=[0;Y_switch];
    Y_switch=abs(Y_switch);

    Y_pswitch(i,:)=Y_switch;
% for n=1:2:(max(HRrequest)-1)
% 
% 
% 
%     Y_pswitch(i,n)=mean(Y_switch([find(HRrequest==n);find(HRrequest==n+1)]));
% end
    
    end
end
end


%get block length for blocks



function BL=getBlockLength(block)

BL=[];
for i=1:length(block)
    if ~isempty(block)
        BL(i)=length(block{i}.LRchoice);
    end
end
end

function optimality=getOptimality(block)


RewardPerPress=[];
for i=1:length(block)
        RewardPerPress(i)=(sum(block{i}.HRreward)+sum(block{i}.LRreward))/(sum(block{i}.HRpress)+sum(block{i}.LRpress));
        optimality(i)=RewardPerPress(i)/(69/200);
end

    
end
