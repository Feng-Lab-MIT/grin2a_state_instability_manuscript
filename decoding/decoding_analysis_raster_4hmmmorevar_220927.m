%20220406 hmm_v2: select blocklength and use only cont1
%20220407 hmm_morevar_2: add more previous choices and previous press N (fix the
%discrete issues)
%20220418 <= smoothing of Neuroninfo  
%20220428 merge HMM 1 & 3 state,
%(get_firing_regressor_per_session.....sm2.m line346=> interested event =16
%(hmm state only)
%20220524 loss the selection criteria => 90% trials should have at
%least one spikes
%20220927 <= include pressn-1 pressn+1

%use only HR trial for reward

%use stricter criteria for non firing blocks (100% trials should have at
%least one spikes

%dividing into only 2 categories

%add next_choice for reward

%use only RS for PFC neurons
%%  1.  Create strings listing where the toolbox and the tutoral data directories are

toolbox_directory_name = 'D:\20211206 Decoding Toolbox\ndt.1.0.4\'  % put name of path to the Neural Decoding Toolbox


%%  3.  Add the toolbox to Matlab's path       

addpath(toolbox_directory_name) 
add_ndt_paths_and_init_rand_generator

%% 4. Use Binned data directly 
%binsize 150ms => 300ms
%steps 50ms??? 
%binned_data (1x N_neurons, no z-score)
%binned_labels (struct 1x1 each one events) 

Neurontype='PL';
Binning=0.001;
fieldname='Initiation';
Smoothing=20;

%binned_data_file_name=strcat('D:\20211206 Decoding Toolbox\',Neurontype,fieldname,'.mat');

timeinibeg=-5;
timeiniend=2;
timerewbeg=-2;
timerewend=7;
edgesini=timeinibeg:Binning:timeiniend;
edgesrewards=timerewbeg:Binning:timerewend;
%%
S=dir('Z:\LeverPressing\ephys');
S2=dir('Z:\Ray data resort\');
S=cat(1,S,S2);

clearvars S2

Neuroninfoall=[];% 
%Regressorall=[];
%%
sessi=1;
for si=1:numel(S)
    
        if ~isempty(dir(strcat('Z:\LeverPressing\ephys\',S(si).name,'\','UnalignedData.mat')))||~isempty(dir(strcat('Z:\Ray data resort\',S(si).name,'\','UnalignedData.mat')))
            if ~isempty(dir(strcat('Z:\LeverPressing\ephys\',S(si).name,'\','RasterData.mat')))
                RasterDatafilepath=strcat('Z:\LeverPressing\ephys\',S(si).name,'\','RasterData.mat');
                UnalignedDatafilepath=strcat('Z:\LeverPressing\ephys\',S(si).name,'\','UnalignedData.mat');
            else
                RasterDatafilepath=strcat('Z:\Ray data resort\',S(si).name,'\','RasterData.mat');
                UnalignedDatafilepath=strcat('Z:\Ray data resort\',S(si).name,'\','UnalignedData.mat');
            end
            
            %for session 11-23, flip the labels
            if ((((strcmp(S(si).name,'2021-11-23_18-56-13')==1||(strcmp(S(si).name,'2022-01-13_23-29-30')==1||strcmp(S(si).name,'2022-01-05_17-23-12')==1))||(strcmp(S(si).name,'2022-01-08_18-47-35')==1))||(strcmp(S(si).name,'2021-12-30_16-55-18')==1))||(strcmp(S(si).name,'2022-02-28_18-30-41')==1))||(strcmp(S(si).name,'2022-03-18_19-44-57')==1)
            
                if strcmp(Neurontype,'MD')==1  
                    [regressor,Neuroninfo,interestedblocktrial] = get_firing_regressor_per_session_initiation_hmm_morevar_220927(RasterDatafilepath,UnalignedDatafilepath,'PL',Binning,fieldname,edgesini,edgesrewards,Smoothing);
                else
                    [regressor,Neuroninfo,interestedblocktrial] = get_firing_regressor_per_session_initiation_hmm_morevar_220927(RasterDatafilepath,UnalignedDatafilepath,'MD',Binning,fieldname,edgesini,edgesrewards,Smoothing);
                end
            else 
                [regressor,Neuroninfo,interestedblocktrial] = get_firing_regressor_per_session_initiation_hmm_morevar_220927(RasterDatafilepath,UnalignedDatafilepath,Neurontype,Binning,fieldname,edgesini,edgesrewards,Smoothing);
            end
            
            if isempty(Neuroninfoall)
                Neuroninfoall=Neuroninfo;
                %Regressorall=regressor;
            else
                Neuroninfoall=cat(2,Neuroninfoall,Neuroninfo);
                %Regressorall=cat(1,Regressorall,regressor);
            end
            

            clearvars Neuroninfo RasterDatafilepath UnalignedDatafilepath regressor

            sessi=sessi+1;
            
        end

end

clearvars si sessi RasterDatafilepath UnalignedDatafilepath


%%
%<= DO NOT resample the trials

%select neurons

%1. select firing rate 491=>390
SelectNeuron_FR=[Neuroninfoall(1:numel(Neuroninfoall)).inverseISIFRcheck]';

%2. select session => 491=> 190
SessTrialCutoff=30; %select non-zero trial session (maybe increase this?)

if strcmp(fieldname,'Initiation')==1
    potentialbinned_data={Neuroninfoall(1:numel(Neuroninfoall)).IniStack};
    SelectNeuron_Sess=zeros(size(SelectNeuron_FR,1),size(SelectNeuron_FR,2));
    for i=1:numel(potentialbinned_data) 
       SelectNeuron_Sess(i,1)=(size(potentialbinned_data{i},1)>SessTrialCutoff)&(sum(sum(potentialbinned_data{i},2)>0)>=0.9*size(potentialbinned_data{i},1)); %select non-zero session with 0.9
    
    end
else
    potentialbinned_data={Neuroninfoall(1:numel(Neuroninfoall)).RewStack};
    SelectNeuron_Sess=zeros(size(SelectNeuron_FR,1),size(SelectNeuron_FR,2));
    for i=1:numel(potentialbinned_data) 
       SelectNeuron_Sess(i,1)=(size(potentialbinned_data{i},1)>SessTrialCutoff)&(sum(sum(potentialbinned_data{i},2)>0)>=0.9*size(potentialbinned_data{i},1)); %select non-zero session with 0.9
    
    end
end
%%

%3. make sure it is RS if it is PL
%select RS neurons (MD:491=>491
load('D:\20220921 decoding\FS_RS_4moresession.mat','RSidx','Num_seq');

SelectNeuron_RS=-1*ones(size(SelectNeuron_FR,1),size(SelectNeuron_FR,2));
if strcmp(Neurontype,'PL')==1
    for ni=1:numel(Neuroninfoall)
        for fsi=1:length(Num_seq)
            if strcmp(Neuroninfoall(ni).Sec(1:6),Num_seq{fsi,1}([3:4,6:7,9:10]))==1
                if strcmp(Neuroninfoall(ni).TT(1:2),'TT')==1
                    if (strcmp(Neuroninfoall(ni).TT(3:end),num2str(Num_seq{fsi,2}))==1)&&(Neuroninfoall(ni).UnitNbr==Num_seq{fsi,3})
                        SelectNeuron_RS(ni,1)=RSidx(fsi,1);
                    end
                    
                else
                    if (strcmp(Neuroninfoall(ni).TT,Num_seq{fsi,2})==1)&&(Neuroninfoall(ni).UnitNbr==Num_seq{fsi,3})
                        SelectNeuron_RS(ni,1)=RSidx(fsi,1);
                    end
                end
    
            end
        end
    end
    
else
    SelectNeuron_RS=1*ones(size(SelectNeuron_FR,1),size(SelectNeuron_FR,2));
end


selectedneuron=((SelectNeuron_FR) &  (SelectNeuron_Sess==1))&(SelectNeuron_RS==1);
selectedneuronid=find(selectedneuron==1);

clearvars potentialbinned_data
%%
potentialbinned_data=[];
%2. DO NOT match trial number
if strcmp(fieldname,'Initiation')==1
    potentialbinned_data={Neuroninfoall(selectedneuron).IniStack};
else
    for i=1:length(selectedneuronid)
        potentialbinned_data{i}=Neuroninfoall(selectedneuronid(i)).RewStack(strcmp(Neuroninfoall(selectedneuronid(i)).current_choice,'1')==1,:);
    end
end

%%
fn=fieldnames(Neuroninfoall);

if strcmp(fieldname,'Initiation')==1
    for fni=10:length(fn)
        potentialbinned_labels.(fn{fni})={Neuroninfoall(selectedneuron).(fn{fni})};
    end
else
    for fni=11:length(fn) %for reward: start from next to current choice
        for i=1:length(selectedneuronid)
            potentialbinned_labels.(fn{fni}){i}=Neuroninfoall(selectedneuronid(i)).(fn{fni})(strcmp(Neuroninfoall(selectedneuronid(i)).current_choice,'1')==1);
        end
    end
end

%%
raster_data_directory_name=strcat('D:\20220921 decoding\','PL_Ini_w4moreRaysess2','\');

mkdir PL_Ini_w4moreRaysess2


%%
for ineuron=1:numel(potentialbinned_data)
    raster_data=potentialbinned_data{ineuron};
    
    raster_site_info.session_ID=Neuroninfoall(selectedneuronid(ineuron)).Sec;
    raster_site_info.recording_channel=Neuroninfoall(selectedneuronid(ineuron)).TT;
    raster_site_info.unit=Neuroninfoall(selectedneuronid(ineuron)).UnitNbr;
    
    %labels
    if strcmp(fieldname,'Initiation')==1
        for fni=10:length(fn)
            raster_labels.(fn{fni})=potentialbinned_labels.(fn{fni}){ineuron};
        end
    else
        for fni=11:length(fn)
            raster_labels.(fn{fni})=potentialbinned_labels.(fn{fni}){ineuron};
        end
    end
    
    rasterfilename=strcat(Neurontype,Neuroninfoall(selectedneuronid(ineuron)).Sec,Neuroninfoall(selectedneuronid(ineuron)).TT,'_',num2str(Neuroninfoall(selectedneuronid(ineuron)).UnitNbr),'.mat');
    save([raster_data_directory_name rasterfilename],'raster_data','raster_labels','raster_site_info');

end

clearvars raster_data raster_labels rasterfilename raster_site_info potentialbinned_data potentialbinned_labels

%% data bias <<= around 6:4 

%load('D:\20211206 Decoding Toolbox\Rasterdata\MD210901TT24_3.mat', 'raster_labels');
%sum(strcmp(raster_labels.comingpress,'1')==1)
%sum(strcmp(raster_labels.comingpress,'-1')==1)


%%  4.  Bin the data

save_prefix_name = strcat('Binned_',Neurontype,fieldname,'4more_sess2');
bin_width = 600; 
step_size = 200;  

binned_data_file_name = create_binned_data_from_raster_data(raster_data_directory_name, save_prefix_name, bin_width, step_size);
 

%%
% binned_data=potentialbinned_data;
% 
% binned_labels.comingpress=potentialbinned_labels.comingpress;
% 
% binned_site_info=[];
% 
% % save extra information about the bin_width, sampling_interval, etc.
% %binned_site_info.binning_parameters.raster_file_directory_name = raster_file_directory_name;
% binned_site_info.binning_parameters.bin_width = Binning*1000;
% binned_site_info.binning_parameters.sampling_interval = Binning*1000; %sampling_interval;%step sizes
% 
% % binned_site_info.binning_parameters.start_time = start_time;
% % binned_site_info.binning_parameters.end_time = end_time;
% % 
% % binned_site_info.binning_parameters.the_bin_start_times = the_bin_start_times;
% % binned_site_info.binning_parameters.the_bin_widths = the_bin_widths;
% 
% save(binned_data_file_name,'binned_data','binned_labels','binned_site_info');


%%  5.  Calculate how many times each stimulus has been shown to each neuron <<??????
%binned_data_file_name='Binned_PLInitiationmorevarRS_600ms_bins_200ms_sampled.mat';
load(binned_data_file_name);  % load the binned data

fn2=fieldnames(binned_labels);

k=1;
for fni=1:length(fn2)
    for i = 0:60
        inds_of_sites_with_at_least_k_repeats = find_sites_with_k_label_repetitions(binned_labels.(fn2{fni}), i);
        num_sites_with_k_repeats(k,i + 1) = length(inds_of_sites_with_at_least_k_repeats);
    end
    k=k+1;
end

%20220406for hmm, MD => 26 are 187/187
%20220406for hmm, PL => 26 are 168/168

%20220406for hmm2, PL => 13 are 128/147 (new data:113, all data:147
%neurons)

%20220406for hmm2, MD => 13 are 127/147 (new data:107, all data:147
%neurons)


%%  Begin the decoding analysis  %%

% 
% 
% for fni=1
% %%  6.  Create a datasource object
% 
% % we will use object identity labels to decode which object was shown (disregarding the position of the object)
% %specific_binned_labels_names = 'past_choice';
% specific_binned_labels_names = fn2{fni};
% 
% % use 20 cross-validation splits (which means that 19 examples of each object are used for training and 1 example of each object is used for testing)
% num_cv_splits = 10; 
% 
% % create the basic datasource object
% ds = basic_DS(binned_data_file_name, specific_binned_labels_names,  num_cv_splits);
% 
% 
% 
% % other useful options:
% 
% % if using the Poison Naive Bayes classifier, load the data as spike counts by setting the load_data_as_spike_counts flag to 1
% %ds = basic_DS(binned_data_file_name, specific_binned_labels_names,  num_cv_splits, 1);
% 
% % can have multiple repetitions of each label in each cross-validation split (which is a faster way to run the code that uses most of the data)
% ds.num_times_to_repeat_each_label_per_cv_split = 1;
% 
%  % optionally can specify particular sites to use
% %ds.sites_to_use = find_sites_with_k_label_repetitions(the_labels_to_use, num_cv_splits);  
% ds.sites_to_use = find_sites_with_k_label_repetitions(binned_labels.(specific_binned_labels_names), num_cv_splits*ds.num_times_to_repeat_each_label_per_cv_split); 
% 
% % can do the decoding on a subset of labels
% %ds.label_names_to_use =  {'kiwi', 'flower', 'guitar', 'hand'};
% 
% %%
% %check ds <<< 
% [XTr_all_time_cv YTr_all XTe_all_time_cv YTe_all] = get_data(ds);
% 
% 
% %%   7.  Create a feature preprocessor object
% 
% % create a feature preprocess that z-score normalizes each feature
% the_feature_preprocessors{1} = zscore_normalize_FP;  
% 
% 
% % other useful options:   
% 
% % can include a feature-selection features preprocessor to only use the top k most selective neurons
% % fp = select_or_exclude_top_k_features_FP;
% % fp.num_features_to_use = 25;   % use only the 25 most selective neurons as determined by a univariate one-way ANOVA
% % the_feature_preprocessors{2} = fp;
% 
% 
% 
% 
% %%  8.  Create a classifier object 
% 
% % select a classifier
% %the_classifier = max_correlation_coefficient_CL;
% 
% 
% % other useful options:   
% 
% % use a poisson naive bayes classifier (note: the data needs to be loaded as spike counts to use this classifier)
% %the_classifier = poisson_naive_bayes_CL;  
% 
% % use a support vector machine (see the documentation for all the optional parameters for this classifier)
% the_classifier = libsvm_CL;
% 
% 
% %%  9.  create the cross-validator 
% 
% 
% the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);  
% 
% the_cross_validator.num_resample_runs = 40;  % usually more than 2 resample runs are used to get more accurate results, but to save time we are using a small number here
% 
% 
% % other useful options:   
% 
% % can greatly speed up the run-time of the analysis by not creating a full TCT matrix (i.e., only trainging and testing the classifier on the same time bin)
% % the_cross_validator.test_only_at_training_times = 1;  
% 
% 
% 
% 
% %%  10.  Run the decoding analysis   
% 
% % if calling the code from a script, one can log the code so that one can recreate the results 
% %log_code_obj = log_code_object;
% %log_code_obj.log_current_file; 
% 
% 
% % run the decoding analysis 
% DECODING_RESULTS = the_cross_validator.run_cv_decoding; 
% 
% 
% 
% %%  11.  Save the results
% 
% % save the results
% save_file_name = strcat(Neurontype,fieldname,specific_binned_labels_names,'hmm20220406');
% save(save_file_name, 'DECODING_RESULTS'); 
% 
% % if logged the code that was run using a log_code_object, save the code
% %LOGGED_CODE = log_code_obj.return_logged_code_structure;
% %save(save_file_name, '-v7.3', 'DECODING_RESULTS', 'LOGGED_CODE'); 
% 
% 
% 
% %%  12.  Plot the basic results
% 
% 
% % which results should be plotted (only have one result to plot here)
% result_names{1} = save_file_name;
% 
% % create an object to plot the results
% plot_obj = plot_standard_results_object(result_names);
% 
% % put a line at the time when the stimulus was shown
% if strcmp(fieldname,'Initiation')
%     plot_obj.significant_event_times = -timeinibeg*1000;
% else
%     plot_obj.significant_event_times = -timerewbeg*1000;
% end
% plot_obj.errorbar_type_to_plot=1;
% plot_obj.errorbar_file_names=result_names;
% 
% % optional argument, can plot different types of results
% %plot_obj.result_type_to_plot = 3;  % for example, setting this to 2 plots the normalized rank results
% 
% 
% plot_obj.plot_results;   % actually plot the results
% ylim([30 90]);
% %saveas(gcf,strcat(Neurontype,'_RS_',fieldname(1:3),specific_binned_labels_names,'sp10r40'),'fig');
% 
% %keyboard;
% close all
% end
% %% fix the_time_interval
% 
% 
% 
% 
% 
% 
% %%  13.  Plot the TCT matrix
% figure()
% plot_obj = plot_standard_results_TCT_object(save_file_name);
% 
% plot_obj.significant_event_times = 0;   % the time when the stimulus was shown
% 
% 
% % optional parameters when displaying the TCT movie
% %plot_obj.movie_time_period_titles.title_start_times = [-500 0];
% %plot_obj.movie_time_period_titles.title_names = {'Fixation Period', 'Stimulus Period'}
% 
% plot_obj.plot_results;  % plot the TCT matrix and a movie showing if information is coded by a dynamic population code
% 

