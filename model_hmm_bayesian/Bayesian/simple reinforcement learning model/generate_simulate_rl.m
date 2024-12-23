function [actions,pressn,state]=generate_simulate_rl(alpha,hr_m0,nblock)



    pressntest=[1,1,1,1:1:70];
    pressnseq=[1:100];
    
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
    %generate actions
    
    phr_hr_per_trial=phr_large_all;
    
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
                    
                else
                    phr=phr_hr_per_trial(statei+4);

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

