function [NLL]=fitBayesianModelextreme_alpha_noinertia_fromPHR_231130(actions,pressn,phr_hr_per_trial)

    
    %%
    %calculate NLL
    
    %phr_hr_per_trial=mean(phr_large_all,2);
    
    NLL=0;
    for block=1:length(actions)
        
        loglikelihood=zeros(length(actions{block}),1);
            
        numberofhr1s=1;
        phr_hr_pre_hr=phr_hr_per_trial(1);
            for trial=1:length(actions{block}) %number of trials
                
                %determine inertia term based on pressn


                if actions{block}(trial)>0
                    if pressn{block}(trial)==1
                        loglikelihood(trial)=log(phr_hr_per_trial(min([numberofhr1s,4]),1));
                        phr_hr_pre_hr=phr_hr_per_trial(min([numberofhr1s,4]),1);
                        numberofhr1s=numberofhr1s+1;
                    else
                        loglikelihood(trial)=log(phr_hr_per_trial(pressn{block}(trial)+3));
                        phr_hr_pre_hr=phr_hr_per_trial(pressn{block}(trial)+3);
                    end
                    
                else

                    loglikelihood(trial)=log(1-phr_hr_pre_hr);

                end
                %fprintf('Trial %d\n',trial);

            end
            phr_hr_pre_hr=[];

        NLL=NLL-sum(loglikelihood,'omitnan');

    end
    


end

