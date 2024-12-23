function [phr_hr_per_trial,phr_large_all]=generate_phr_BayesianModelextreme_alpha_noinertia(alpha,noisyfactorhr,hr_m0,hr_std0,lr_std)

%20220311 simulate model => get phr for each condition first
%consider actions and pressn from block start

%20220315 get p(hr large) from pure simulation <- don't worry about tpdf or
%integration

%20220316 remove the constrains of hrreward/cost>0 

%20220406 remove lr_m 

%20220628: remove inertia if hr>inertia thres

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
    phr_large_all=zeros(73,200);

    hr_b_all=zeros(73,200);
    hr_m_all=zeros(73,200);
    hr_std_all=zeros(73,200);

    for resample=1:200 %number of trials

        %hr_a=hr_a0;
        %hr_b=hr_b0;
        %hr_k=hr_k0;
        hr_m=hr_m0;
        hr_std=hr_std0;

        for trial=1:length(pressntest)


            %likelihood of choosing action
            %phr_large=0;

            %this is slower but more accurate
            %for xidx=1:length(x) %too slow using trapz, using results from gaussian intergral instead
            %    phr_large=phr_large+0.1*pdf('tLocationScale',x(xidx),hr_m, hr_std,2*hr_a).*cdf('tLocationScale',x(xidx),3*lr_m,3*lr_std,2*lr_a);
            %end


            %issue fixed! (add 1/hr_std to pdf)
%             for xidx=1:length(x) 
%                 phr_large=phr_large+0.5*(1/hr_std)*tpdf((x(xidx)-hr_m)/hr_std,2*hr_a).*tcdf((x(xidx)-3*lr_m)/(3*lr_std),2*lr_a); %issue solved!
%             end

%             phr_large_all(trial,resample)=1-phr_large;
            
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


    
    phr_hr_per_trial=mean(phr_large_all,2);
    


end

