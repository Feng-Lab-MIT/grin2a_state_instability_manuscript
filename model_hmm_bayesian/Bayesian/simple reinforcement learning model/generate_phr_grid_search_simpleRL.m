%version 2 20231128, reduce the parameter


alphalist=[0.1:0.9];


%%

NLLall=zeros(9,1);
for m=1

   alpha=alphalist(m);
   hr_m0=1;

    
   [NLL]=fitModel_RL(alpha,hr_m0,actions,pressn)
   
   NLLall(m,1)=NLL;

    
end


