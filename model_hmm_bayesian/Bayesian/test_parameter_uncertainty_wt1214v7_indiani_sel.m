%bootstrap cells
%20231206 v2: try sample 200 blocks per round instead of 100 
%v3: try 300 blocks 
%v5: use more alpha 
%v6: use even more alpha and use 200 blocks
%v7: use different data set
%1214v7: use more noisyfactor, more hr std.

%20231206 remove : % RParray=cell2mat(mat_block);% pressn={RParray.HRpress}; % actions={RParray.HRreward};
%%

%run this code "combine_v6_v7_grid_search.m" to get PHR_Large_mean_combined_sub;
load("Z:\20231027 lever pressing paper figure\Bayesian\20231101 bayesian model\phr_grid_search_combined_sub_v7_20231214.mat","PHR_Large_mean_combined");

PHR_Large_mean_combined=PHR_Large_mean_combined([3,5:8,10,11,14:15],[1,3:9],[2:9],[2:6],:,:);

%%

%load('NLLcon_finmin_wt_v5_20230422.mat', 'WT_sub_allcells');
%load('NLLcon_finmin_grin2a_v5_20230423.mat', 'grin2a_sub_allcells');
%load('D:\20231027 lever pressing paper figure\Bayesian\20231101 bayesian model\Sum_Data_Seperate_Animals_allcellsflat_231208_noselect.mat','grin2a_sub','WT_sub');
load('Z:\20231027 lever pressing paper figure\Bayesian\20231101 bayesian model\Sum_Data_Seperate_Animals_allcellsflat_231208.mat','grin2a_sub','WT_sub');
%load('PHR_grin_search_coarse.mat')

mat_sub=WT_sub;

bestparameter=zeros(numel(fieldnames(mat_sub)),4);

aniname=fieldnames(mat_sub);

for fi=1:numel(fieldnames(mat_sub))

    mat_block=mat_sub.(aniname{fi});
    
    for blocki=1:length(mat_block)
        if sum(strcmp(fieldnames(mat_block{blocki}),'blockOnset'))>0
            if mat_block{blocki}.blockOnset>0
                actions{blocki}=mat_block{blocki}.HRreward(mat_block{blocki}.blockOnset:end)>0;
                pressn{blocki}=mat_block{blocki}.HRpress(mat_block{blocki}.blockOnset:end);
            elseif mat_block{blocki}.blockOnset==0
                blockonset=find(mat_block{blocki}.HRreward(1:find(mat_block{blocki}.HRpress==2)-4)>0,1,'last');
                actions{blocki}=mat_block{blocki}.HRreward(blockonset:end)>0;
                pressn{blocki}=mat_block{blocki}.HRpress(blockonset:end);
                
            end
        else
            blockonset=find(mat_block{blocki}.HRreward(1:find(mat_block{blocki}.HRpress==2)-4)>0,1,'last');
            actions{blocki}=mat_block{blocki}.HRreward(blockonset:end)>0;
            pressn{blocki}=mat_block{blocki}.HRpress(blockonset:end);
        end
    end
    
    
    
    
    %% sample with replacement 
    NLL=zeros(9,8,8,5);

     for k=1:8 %hr std 0
       for l=1:5
          for i = 1:9
             for j = 1:8%min(floor(k*2.5),7) %noisyfactor
                 NLL(i,j,k,l)=fitBayesianModelextreme_alpha_noinertia_fromPHR_231130(actions,pressn,reshape(PHR_Large_mean_combined(i,j,k,l,1,:),[],1));
             end
          end
       end
     end
     
    
     ind=find(-NLL==max(-NLL(:)));
     [a,b,c,d]=ind2sub(size(NLL),ind);
     bestparameter(fi,:)=[a,b,c,d];
     

end

%% save

save('individual_animal_wt_noselect_231214sel.mat','bestparameter','WT_sub');