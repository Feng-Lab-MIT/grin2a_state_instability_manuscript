function [NLL]=fitModel_RL(alpha,hr_m0,actions,pressn)




    pressntest=[1,1,1,1:1:70];
    phr_large_all=zeros(73,1);

    lr_m=6;
    hr_m=hr_m0;

    for trial=1:length(pressntest)
        

        hrvaluecost=3/hr_m;
        lrvaluecost=1/lr_m;
        
        if (hrvaluecost>lrvaluecost)
            phr_large_all(trial)=1-10^(-100);
        elseif (hrvaluecost<lrvaluecost)
            phr_large_all(trial)=10^(-100);
        else
            phr_large_all(trial)=0.5;
        end


        %================
        %how to code noisy <<<<<
        %================
        %if experienced_hrrequest+2*(random.random()-0.5)*self.noisyfactorhr*experienced_hrrequest>0:
        %noisyhr=experienced_hrrequest+2*(random.random()-0.5)*self.noisyfactorhr*experienced_hrrequest

        experienced_hrrequest=pressntest(trial);

        %experienced_hrrequest=max([0,experienced_hrrequest+2*(rand(1)-0.5)*noisyfactorhr*experienced_hrrequest]);
        

        %================
        %update HRrequest_belief distribution
        %================


        hr_m_new=hr_m+alpha*(experienced_hrrequest-hr_m);

        hr_m=hr_m_new;


    end



    
    %%
    %calculate NLL
    
    phr_hr_per_trial=phr_large_all;
    
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
                    if (1-phr_hr_pre_hr)==0
                        loglikelihood(trial)=log(10^(-100));
                    else
                        loglikelihood(trial)=log(1-phr_hr_pre_hr);
                    end

                end
                %fprintf('Trial %d\n',trial);

            end
            phr_hr_pre_hr=[];

        NLL=NLL-sum(loglikelihood,'omitnan');

    end
    


end

