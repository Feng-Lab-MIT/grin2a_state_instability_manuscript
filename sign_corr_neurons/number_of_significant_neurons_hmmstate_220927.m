%v2: 20220203 include 22-01-05 session and 22-01-13 2 sessions (including
%changing the FS/RS labeling <- reclustering)
%v3: 20220204: use corr-coef instead of regression (code changes in
%regression_with_model_parameters_v6_5_selblock_ambbehunc_nofig.m
%20220314: add 2 more sections and use new version of regression_with_model_parameters_v6_5_selblock_ambbehunc_nofig.m
%new version fix the HR request to previous HR request instead of upcoming

%20220425: add more previous trials and HMM states

Neurontype='MD';
Neuronsubtype='RS';
Binning=0.5; %window:2 increment:0.5
fieldname='Initiation';

% if strcmp(fieldname,'Initiation')==1
%     interestevent=[1,3,4,5,6,8,15];%4:upcoming HR 5:past HR
% else
%     interestevent=[1,2,4,6,8,15];%upcoming choice
% end

%1: HR/LRchoice (current)
%2: next choice -reward
%3: past choice -ini
%4: t-2 choice
%5: t-3 choice
%6: t-4 choice
%7: t-5 choice
%8: t-6 choice
%9: press number (current) -reward
%10: press number (next) - rew
%11: press number (past) - ini
%12: t-2 hr press number
%13: t-3 hr press number
%14: t-4 hr press number
%15: t-5 hr press number
%16: reward/cost (current) -reward
%17: reward/cost (next) -rew
%18: reward/cost (past) - ini
%19: t-2 reward/cost
%20: t-3 reward/cost
%21: t-4 reward/cost
%22: t-5 reward/cost
%23: hmm states 
%24: behavior uncertainty

timeinibeg=-5;
timeiniend=2;
timerewbeg=-2;
timerewend=7;
Windowsize=2;
Increment=0.5;

edgesini=timeinibeg:Binning:timeiniend;
edgesrewards=timerewbeg:Binning:timerewend;
%%
S=dir('X:\Tingting\LeverPressing\ephys');
S2=dir('X:\Yi_Yun\Ray data resort\');
S=cat(1,S,S2);

clearvars S2

Neuroninfoall=[];% 
RegRsqrIniall=[];
RegRsqrIniSignNeuall=[];
RegRsqrRewall=[];
RegRsqrRewSignNeuall=[];

%Regressorall=[];

%%
sessi=1;
for si=1:numel(S)-3
    
        if ~isempty(dir(strcat('X:\Tingting\LeverPressing\ephys\',S(si).name,'\','UnalignedData.mat')))||~isempty(dir(strcat('X:\Yi_Yun\Ray data resort\',S(si).name,'\','UnalignedData.mat')))
            if ~isempty(dir(strcat('X:\Tingting\LeverPressing\ephys\',S(si).name,'\','RasterData.mat')))
                RasterDatafilepath=strcat('X:\Tingting\LeverPressing\ephys\',S(si).name,'\','RasterData.mat');
                UnalignedDatafilepath=strcat('X:\Tingting\LeverPressing\ephys\',S(si).name,'\','UnalignedData.mat');
            else
                RasterDatafilepath=strcat('X:\Yi_Yun\Ray data resort\',S(si).name,'\','RasterData.mat');
                UnalignedDatafilepath=strcat('X:\Yi_Yun\Ray data resort\',S(si).name,'\','UnalignedData.mat');
            end

            [Neuronfiringrate,RegressRsqr_Ini,RegressRsqr_Ini_SignNeuron,targetini,RegressRsqr_Reward,RegressRsqr_Reward_SignNeuron,targetreward] = regression_with_hmmstate(RasterDatafilepath,UnalignedDatafilepath,Neurontype,Windowsize,Increment,timeinibeg,timeiniend,timerewbeg,timerewend);    

            
            if isempty(Neuroninfoall)
                Neuroninfoall=Neuronfiringrate;
                RegRsqrIniall=RegressRsqr_Ini;
                RegRsqrIniSignNeuall=RegressRsqr_Ini_SignNeuron;
                RegRsqrRewall=RegressRsqr_Reward;
                RegRsqrRewSignNeuall=RegressRsqr_Reward_SignNeuron;
            else
                Neuroninfoall=cat(2,Neuroninfoall,Neuronfiringrate);
                RegRsqrIniall=cat(1,RegRsqrIniall,RegressRsqr_Ini);
                RegRsqrIniSignNeuall=cat(1,RegRsqrIniSignNeuall,RegressRsqr_Ini_SignNeuron);
                RegRsqrRewall=cat(1,RegRsqrRewall,RegressRsqr_Reward);
                RegRsqrRewSignNeuall=cat(1,RegRsqrRewSignNeuall,RegressRsqr_Reward_SignNeuron);
            end
            

            clearvars Neuronfiringrate RasterDatafilepath UnalignedDatafilepath RegressRsqr_Ini RegressRsqr_Ini_SignNeuron RegressRsqr_Reward RegressRsqr_Reward_SignNeuron

            sessi=sessi+1;
            
        end

end

clearvars si sessi RasterDatafilepath UnalignedDatafilepath

%%
%<= DO NOT resample the trials

%1. select neurons

SessTrialCutoff=100; %select non-zero trial session (maybe increase this?)

SelectNeuron_FR=[Neuroninfoall(1:numel(Neuroninfoall)).inverseISIFRcheck]';

%initiation
potentialbinned_data={Neuroninfoall(1:numel(Neuroninfoall)).IniSelData};
SelectNeuron_SessIni=zeros(size(SelectNeuron_FR,1),size(SelectNeuron_FR,2));
for i=1:numel(potentialbinned_data) 
   SelectNeuron_SessIni(i,1)=(size(potentialbinned_data{i},1)>SessTrialCutoff)&(sum(sum(potentialbinned_data{i},2)>0)>=1*size(potentialbinned_data{i},1)); %select non-zero session

end

%reward
potentialbinned_data={Neuroninfoall(1:numel(Neuroninfoall)).RewSelData};
SelectNeuron_SessRew=zeros(size(SelectNeuron_FR,1),size(SelectNeuron_FR,2));
for i=1:numel(potentialbinned_data) 
   SelectNeuron_SessRew(i,1)=(size(potentialbinned_data{i},1)>SessTrialCutoff)&(sum(sum(potentialbinned_data{i},2)>0)>=1*size(potentialbinned_data{i},1)); %select non-zero session

end

SelectNeuron_Sess=SelectNeuron_SessIni&SelectNeuron_SessRew;
%%
%load('D:\20211206 Decoding Toolbox\20211217 FS RS\20220203 add 3 more sessions\Waveform_para_all20220203.mat','RSidx','Num_seq');
%load('D:\20211206 Decoding Toolbox\20211217 FS RS\20220204 add 2 more sessions\Waveform_para_all20220204.mat','RSidx','Num_seq');
load('D:\20220921 decoding\FS_RS_4moresession.mat','RSidx','Num_seq');

%2. make sure it is RS if it is PL 
SelectNeuron_RS=-1*ones(size(SelectNeuron_FR,1),size(SelectNeuron_FR,2));
SelectNeuron_FS=-1*ones(size(SelectNeuron_FR,1),size(SelectNeuron_FR,2));
if strcmp(Neurontype,'PL')==1
    for ni=1:numel(Neuroninfoall)
        for fsi=1:length(Num_seq)
            if strcmp(Neuroninfoall(ni).Sec(1:6),Num_seq{fsi,1}([3:4,6:7,9:10]))==1
                if strcmp(Neuroninfoall(ni).TT(1:2),'TT')==1
                    if (strcmp(Neuroninfoall(ni).TT(3:end),num2str(Num_seq{fsi,2}))==1)&&(Neuroninfoall(ni).UnitNbr==Num_seq{fsi,3})
                        SelectNeuron_RS(ni,1)=RSidx(fsi,1);
                        SelectNeuron_FS(ni,1)=~RSidx(fsi,1);
                    end
                    
                else
                    if (strcmp(Neuroninfoall(ni).TT,Num_seq{fsi,2})==1)&&(Neuroninfoall(ni).UnitNbr==Num_seq{fsi,3})
                        SelectNeuron_RS(ni,1)=RSidx(fsi,1);
                        SelectNeuron_FS(ni,1)=~RSidx(fsi,1);
                    end
                end
    
            end
        end
    end
    
else
    SelectNeuron_RS=1*ones(size(SelectNeuron_FR,1),size(SelectNeuron_FR,2));
end

if strcmp(Neuronsubtype,'FS')~=1
    selectedneuron=((SelectNeuron_FR) &  (SelectNeuron_Sess==1))&(SelectNeuron_RS==1);
    selectedneuronid=find(selectedneuron==1);
else
    selectedneuron=((SelectNeuron_FR) &  (SelectNeuron_Sess==1))&(SelectNeuron_FS==1);
    selectedneuronid=find(selectedneuron==1);
end

clearvars potentialbinned_data

%%
save(strcat(Neurontype,'_corr_coef_neuron_20220927.mat'));

