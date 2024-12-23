%20231201 exapnd from version 2, add some more parameters
%20231206 expand alpha 
%4exp: add -0.5 -1.5 in alpha
%5exp: add -8 in alpha, change hr std0list 1,3,2,4 => 0.4,0.6,0.7,0.8
% not finish running
%5_2: alpha 12:14



%alphalist=[-6,-5,-4,-3,-2,-1,0,0.2,0.4,0.6,0.7,0.8,0.9];
%alphalist=[-2,-1,0,0.2,0.4,0.6,0.8,-3,-4]; %(8,9)
%alphalist=[-2,-1,0,0.2,0.4,0.6,0.8,-3,-4,-5,-6,0.9,-0.5,-1.5,-0.2,-0.4,-0.6,-0.8,-1.2,-1.4,-1.6,-1.8]; %(10,11,12)
alphalist=[-2,-1,0,0.2,0.4,0.6,0.8,-3,-4,-0.5,-1.5,-5,-6,0.9,-8]; %(13)
%noisyfactorhrlist=[0.5:0.1:0.9,0.92,0.94,0.96,0.98,1,2,3,5];
noisyfactorhrlist=[0.5,0.9,1,2,3,4]; %6
%hr_std0list=[0.05 0.1 0.2 0.4,0.5,1,2,3,4,5,6];
%hr_std0list=[0.05 0.1 0.2 0.3,0.5,1,3,2,4]; %(8,9)
hr_std0list=[0.05 0.1 0.2 0.3,0.5,0.4,0.6,0.7,0.8]; %(8,9)

lr_stdlist=[0.5,1,2,3,4,5]; %(6)
hr_m0list=[1];



%%
%parameter alpha=-5 noisyfactorhr=8.000000e-01 hr_m0=0 hr_std0=4 lr_std=5

load('PHR_grin_search_coarse_exp_all.mat','PHR_Large_mean','PHR_Large_all');
%%

for m=1
    for k=1:5 %1:9 %1:7 %8:9 %1:9 %hr std 0
       for l=1:6 %1:5 %lr std 
          for i = 12:14 %15 %1:7 %8:9 %alpha
             for j = 1:6 %6 %1:6 %noisyfactor
            
                    alpha=alphalist(i);
                    noisyfactorhr=noisyfactorhrlist(j);
                    
                    hr_std0=hr_std0list(k);
                    lr_std=lr_stdlist(l);
    
                    hr_m0=hr_m0list(m);
                    
                    [phr_hr_per_trial,phr_large_all]=generate_phr_BayesianModelextreme_alpha_noinertia(alpha,noisyfactorhr,hr_m0,hr_std0,lr_std);
            
                    PHR_Large_mean(i,j,k,l,m,:)=phr_hr_per_trial;
                    PHR_Large_all{i,j,k,l,m}=phr_large_all;
            
                    fprintf('PHR_grin_search alpha=%d noisyfactorhr=%d hr_m0=%d hr_std0=%d lr_std=%d\n',alpha,noisyfactorhr,hr_m0,hr_std0,lr_std);

                end
                save('PHR_grin_search_coarse_exp_alpha12_14_231214.mat')

            end
       end
       

    end
    
end


