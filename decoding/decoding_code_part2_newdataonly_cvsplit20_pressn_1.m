% should have run part1 code first and saved binned data


%%  1.  Create strings listing where the toolbox and the tutoral data directories are
toolbox_directory_name = 'D:\20211206 Decoding Toolbox\ndt.1.0.4\'  % put name of path to the Neural Decoding Toolbox

%%  3.  Add the toolbox to Matlab's path       
addpath(toolbox_directory_name) 
add_ndt_paths_and_init_rand_generator

%% 4. Use Binned data directly

for neurontypei=1:2
    
    if neurontypei==1
        Neurontype='MD';
        Binning=0.001;
        fieldname='Initiation';
        binned_data_file_name='Binned_MDInitiation4more_sess2_600ms_bins_200ms_sampled.mat';

    else
        Neurontype='PL';
        Binning=0.001;
        fieldname='Initiation';
        binned_data_file_name='Binned_PLInitiation4more_sess2_600ms_bins_200ms_sampled.mat';

    end
    
    %1: correct or not (1 or 0)
    %2: poke port (4 possibility)
    %3: context (1 or 0)
    %4: left or right
    %5: ifswitch (1 or 0) if at switch states


    %%  5.  Calculate how many times each stimulus has been shown to each neuron <<??????

    Binorigin=load(strcat('D:\20220921 decoding\',binned_data_file_name));  % load the binned data

    fn2=fieldnames(Binorigin.binned_labels);


    for fni=1:length(fn2)
        inds_of_sites_with_at_least_k_repeats = find_sites_with_k_label_repetitions(Binorigin.binned_labels.(fn2{fni}), 20); %find site with at least 10 repeats
        ind_of_sites_to_use{fni}=inds_of_sites_with_at_least_k_repeats;
        num_of_available_site(fni)=length(inds_of_sites_with_at_least_k_repeats);
    end
    
    %use ind_of_sites_to_use for the minimum labels
    ind_of_sites_to_use_final=ind_of_sites_to_use{find(num_of_available_site==min(num_of_available_site),1,'first')};


    %% choose some number of neurons

    % #of run per condi x time frames x # of labels 
    mean_accuracy_80=zeros(50,size(Binorigin.binned_data{1},2),length(fn2));

    %===========================
    % number of neurons to use << change here 
    %===========================
    number_neuron_to_use=120; %check minimum of length(ind_of_sites_to_use)

        
        for validruni=1:20 %intended number trial-1

            %20220518: select from sites with at least 5 repeats
            %selectedneurons_seq=randperm(numel(Binorigin.binned_data),number_neuron_to_use); %20220209 should use randperm instead to avoid repeating the same neurons
            selectedneurons_seq=ind_of_sites_to_use_final(randperm(length(ind_of_sites_to_use_final),number_neuron_to_use));

            
            %get binned_data
            binned_data2={Binorigin.binned_data{[selectedneurons_seq]}};

            %get binned labels
            labelfield=fieldnames(Binorigin.binned_labels);
            for labeli=1:length(labelfield)
                binned_labels2.(labelfield{labeli})={Binorigin.binned_labels.(labelfield{labeli}){[selectedneurons_seq]}};
            end

            %get binned site infor
            binned_site_info2.session_ID={Binorigin.binned_site_info.session_ID{[selectedneurons_seq]}};
            binned_site_info2.recording_channel=Binorigin.binned_site_info.recording_channel([selectedneurons_seq]);
            binned_site_info2.unit=Binorigin.binned_site_info.unit([selectedneurons_seq]);
            binned_site_info2.binning_parameters=Binorigin.binned_site_info.binning_parameters; 

            binned_data=binned_data2;
            binned_labels=binned_labels2;
            binned_site_info=binned_site_info2;

            
            binned_data_file_name2=strcat(binned_data_file_name(1:end-11),'N',num2str(number_neuron_to_use),'run',num2str(validruni),'subsampled_cv20.mat');
            %binned_data_file_name2=strcat('Binned_MDRewardmorevarRS_600ms_bins_200ms_N',num2str(number_neuron_to_use),'subsampled.mat');

            save(binned_data_file_name2,'binned_data','binned_labels','binned_site_info');


            for fni=1:length(fn2)
            %%  Begin the decoding analysis  %%


                %%  6.  Create a datasource object

                % we will use object identity labels to decode which object was shown (disregarding the position of the object)
                specific_binned_labels_names = fn2{fni};

                % use 20 cross-validation splits (which means that 19 examples of each object are used for training and 1 example of each object is used for testing)
                num_cv_splits = 20; %change to 3 <20220205

                % create the basic datasource object
                ds = basic_DS(binned_data_file_name2, specific_binned_labels_names,  num_cv_splits);

                % other useful options:

                % if using the Poison Naive Bayes classifier, load the data as spike counts by setting the load_data_as_spike_counts flag to 1
                %ds = basic_DS(binned_data_file_name, specific_binned_labels_names,  num_cv_splits, 1);

                % can have multiple repetitions of each label in each cross-validation split (which is a faster way to run the code that uses most of the data)
                %ds.num_times_to_repeat_each_label_per_cv_split = 1;

                 % optionally can specify particular sites to use
                %ds.sites_to_use = find_sites_with_k_label_repetitions(the_labels_to_use, num_cv_splits);  
                ds.sites_to_use = find_sites_with_k_label_repetitions(binned_labels.(specific_binned_labels_names), num_cv_splits*ds.num_times_to_repeat_each_label_per_cv_split); 

                % can do the decoding on a subset of labels
                %ds.label_names_to_use =  {'kiwi', 'flower', 'guitar', 'hand'};

                %%
                %check ds <<< 
                [XTr_all_time_cv YTr_all XTe_all_time_cv YTe_all] = get_data(ds);


                %%   7.  Create a feature preprocessor object

                % create a feature preprocess that z-score normalizes each feature
                the_feature_preprocessors{1} = zscore_normalize_FP;  
                % other useful options:   

                % can include a feature-selection features preprocessor to only use the top k most selective neurons
                % fp = select_or_exclude_top_k_features_FP;
                % fp.num_features_to_use = 25;   % use only the 25 most selective neurons as determined by a univariate one-way ANOVA
                % the_feature_preprocessors{2} = fp;

                %%  8.  Create a classifier object 

                % select a classifier
                %the_classifier = max_correlation_coefficient_CL;

                % other useful options:   

                % use a poisson naive bayes classifier (note: the data needs to be loaded as spike counts to use this classifier)
                %the_classifier = poisson_naive_bayes_CL;  

                % use a support vector machine (see the documentation for all the optional parameters for this classifier)
                the_classifier = libsvm_CL;


                %%  9.  create the cross-validator 

                the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);  

                the_cross_validator.num_resample_runs = 10;  % usually more than 2 resample runs are used to get more accurate results, but to save time we are using a small number here

                % other useful options:   

                % can greatly speed up the run-time of the analysis by not creating a full TCT matrix (i.e., only trainging and testing the classifier on the same time bin)
                the_cross_validator.test_only_at_training_times = 1;  

                %%  10.  Run the decoding analysis   

                % if calling the code from a script, one can log the code so that one can recreate the results 
                %log_code_obj = log_code_object;
                %log_code_obj.log_current_file; 

                % run the decoding analysis 
                DECODING_RESULTS = the_cross_validator.run_cv_decoding; 

                %%  11.  Save the results

                % save the results
                mean_accuracy_100(validruni,:,fni)=DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results';
            end
                
                
        end %validrun

        
    accuracyfilename=strcat(binned_data_file_name(1:end-11),'meanaccuracy_cv20.mat');
    save(accuracyfilename);
    clear;
    
    
end
