function [actions,pressn,state]=generate_simulate(alpha,noisyfactorhr,hr_m0,hr_std0,lr_std,nblock)

%20220311 simulate model => get phr for each condition first
%consider actions and pressn from block start

%20220315 get p(hr large) from pure simulation <- don't worry about tpdf or
%integration

%20220316 remove the constrains of hrreward/cost>0 

%20220406 remove lr_m 

%20220624: for optimization using 7 parameters
%parameters from model_fitting_draft_nogrammextre_inertia_exp.m 
%alpha,noisyfactorhr,hr_m0,hr_std0,lr_0,lr_std

%use prob from simulation
%     alpha0 
%     noisyfactorhr

    %input actions=multiple actions in cell array
    %input pressn=multiple pressn in cell array
    
    %x=-200:0.5:200;

    %hr_a0=1.6355;
    %hr_b0=15.2426;
    %hr_k0=3.2710;
    hr_a=0.5;
    hr_k=1-alpha;
    %hr_std0=3.5;
    %hr_m0=1.5028;
    %hr_std=( hr_b*(hr_k+1)/(hr_a*hr_k) )^(0.5);

    %lr value is the average from averging from "D:\20220113 model revision II\previous
    %version\results from original code _lrdistoutput\BayesianRLmodels_normalgamma_full_noquit_noisy_no_noisealpha_al_dist_logbatchN1.0alphaadj1.0hrnoise1.8lrnoise0.8memlim7.mat'
    %LR constant 6, std 3
    
    lr_a=0.5;
    lr_m=6;%lr_m0;
    %lr_a=1.23;
    %lr_m=6;
    %lr_std=2.8;


    pressntest=[1,1,1,1:1:70];
    pressnseq=[1:100];
    
    phr_large_all=zeros(73,200);

    hr_b_all=zeros(73,200);
    hr_m_all=zeros(73,200);
    hr_std_all=zeros(73,200);

    for resample=1:200 %number of trials

        hr_m=hr_m0;
        hr_std=hr_std0;

        for trial=1:length(pressntest)


%20220315 phr_large from pure simulation    
            pdhr=makedist('tLocationScale','mu',hr_m,'sigma',hr_std,'nu',2*hr_a);
            pdlr=makedist('tLocationScale','mu',lr_m,'sigma',lr_std,'nu',2*lr_a);

            hrvaluecost=3./random(pdhr,[100000,1]);
            lrvaluecost=1./random(pdlr,[100000,1]);
            
            %markdown 202220316
            %hrvaluecost(hrvaluecost<0)=0;
            %lrvaluecost(lrvaluecost<0)=0;
            phr_large_all(trial,resample)=sum((hrvaluecost-lrvaluecost)>0)/100000;

            %update hr_value and hr_velo in hr trial

            %================
            %how to code noisy <<<<<
            %================
            %if experienced_hrrequest+2*(random.random()-0.5)*self.noisyfactorhr*experienced_hrrequest>0:
            %noisyhr=experienced_hrrequest+2*(random.random()-0.5)*self.noisyfactorhr*experienced_hrrequest

            experienced_hrrequest=pressntest(trial);

            experienced_hrrequest=max([0,experienced_hrrequest+2*(rand(1)-0.5)*noisyfactorhr*experienced_hrrequest]);

            %================
            %update HRrequest_belief distribution
            %================

            %k_new=(1-self.alpha)*k+1          
            %a_new=(1-self.alpha)*a+(1/2) 
            %b_new=(1-self.alpha)*b+k*((hrrequest-mu)**2)/(2*(k+1))
            %mu_new= (k_new*mu+hrrequest)/(k_new+1)
            %std_new= ( b_new*(k_new+1)/(a_new*k_new) )**0.5

            %hr_k=1;
            %hr_a=0.5;
            hr_b_new=hr_k*(((experienced_hrrequest-hr_m)^2)/(2*(hr_k+1)));
            hr_m_new=(hr_k*hr_m+experienced_hrrequest)/(hr_k+1);
            hr_std_new=( hr_b_new*(hr_k+1)/(hr_a*hr_k) )^0.5;

            hr_b=hr_b_new;
            hr_m=hr_m_new;
            hr_std=hr_std_new;

            hr_b_all(trial,resample)=hr_b;
            hr_m_all(trial,resample)=hr_m;
            hr_std_all(trial,resample)=hr_std;


        end

    end %end of phr simulation

    
    %%
    %generate actions
    
    phr_hr_per_trial=mean(phr_large_all,2);
    
    NLL=0;
    for block=1:nblock
           
        phr_hr_pre_hr=phr_hr_per_trial(1);

        trial=1;
        statei=0;
        blockend=0;
        numberofconseqhr1s=0;
        numberofconseqlr=0;
        
        actions{block}=zeros(100,1);
        pressn{block}=zeros(100,1);
        state{block}=zeros(100,1);

        pressn{block}(1)=1;
        while blockend==0
            
%             %determine inertia term based on pressn
%             if pressn{block}(trial)>inertiathres
%                 inertia_pr=0;
%             else
%                 inertia_pr=inertia;
%             end
                
            %1. get phr after considering inertia %phr should from previous
            %states

%             if trial==1
%                 phr_inertia=0.5;
%             
%             else

                if statei==0                        
                    phr=phr_hr_per_trial(min([numberofconseqhr1s+1,4]),1);
                    phr_hr_pre_hr=phr;
                    
                else
                    phr=phr_hr_per_trial(statei+4);
                    phr_hr_pre_hr=phr;
                end
                
%             end


            %2. get phr after considering inertia
            x=rand;

            if x<=phr
                actions{block}(trial)=1;                          
            else
                actions{block}(trial)=0; 
            end
            
            %3. get pressn for current state
            state{block}(trial)=statei;
            pressn{block}(trial)=pressnseq(statei+1);
            
            
            %4. update state 
            
            %determine if state move
            if actions{block}(trial)==1
                
                numberofconseqlr=0;
                if statei==0
                    numberofconseqhr1s=numberofconseqhr1s+1;
                    if numberofconseqhr1s>=4
                        statei=1;                        
                    else
                        statei=0;                        
                    end
                else
                    statei=statei+1;
                end
                
            else
                numberofconseqlr=numberofconseqlr+1;
                numberofconseqhr1s=0;
                
                if numberofconseqlr>=6
                    %reset block
                    blockend=1;
                end
   
            end


            %5. determine if terminate
            
            if trial>=100
                blockend=1;
                break
            end
            
            trial=trial+1;
                    

        end%while
        phr_hr_pre_hr=[];
        actions{block}=actions{block}(1:trial-1,:);
        pressn{block}=pressn{block}(1:trial-1,:);
        state{block}=state{block}(1:trial-1,:);


    end
    


end

