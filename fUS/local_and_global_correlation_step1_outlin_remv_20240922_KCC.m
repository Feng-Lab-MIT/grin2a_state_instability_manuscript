%20240922 resample -> KCC

%1. get raw trace from each slices
%2. get selected time window <<<
%3. interpolate => filter  <<<
%4. rotation and align

S=dir(fullfile('\\fenglab03\yiyun\20230515 fUS no dex\local_connectivity_outliner_remove\','*resamplengf1011.mat'));
%S=dir(fullfile('\\fenglab03\yiyun\20240820 fUS\','*resamplengf1017filt.mat')); 
%run ^this first then the previous one

%% KCC
for si=1:numel(S)
    
    load(strcat('\\fenglab03\yiyun\20230515 fUS no dex\local_connectivity_outliner_remove\',S(si).name));
    %load(strcat('\\fenglab03\yiyun\20240820 fUS\',S(si).name)); %run this first then the previous one
    
    KCC_local=zeros(size(datarotate_timecorrect_gf,1),size(datarotate_timecorrect_gf,2),size(datarotate_timecorrect_gf,3));
    Rank=zeros(size(datarotate_timecorrect_gf,1),size(datarotate_timecorrect_gf,2),size(datarotate_timecorrect_gf,3),size(datarotate_timecorrect_gf,4));
    
    % convert datarotate_timecorrect_gf to rank
    


    for xi=1:size(datarotate_timecorrect_gf,1)
        for yi=1:size(datarotate_timecorrect_gf,2)
            for zi=1:size(datarotate_timecorrect_gf,3)

                rawtrace=reshape(datarotate_timecorrect_gf(xi,yi,zi,:),1,[]);
                [~,b]=sort(rawtrace);
                %rank 
                b(2,:)=[1:length(b)];
                r=sortrows(b'); 
                
                Rank(xi,yi,zi,:)=r(:,2)';
                
            end
        end
    end
    
    for xi=1:size(datarotate_timecorrect_gf,1)
        for yi=1:size(datarotate_timecorrect_gf,2)
            for zi=1:size(datarotate_timecorrect_gf,3)   
    

                Nearborcor=[xi,yi,zi]+[
                [0,-1,-1];
                [0,-1,0];
                [0,-1,1];

                [0,0,1];
                [0,0,-1];

                [0,1,-1];
                [0,1,0];
                [0,1,1];

                [1,-1,-1];
                [1,-1,0];
                [1,-1,1];

                [1,0,-1];
                [1,0,0];
                [1,0,1];

                [1,1,-1];
                [1,1,0];
                [1,1,1];

                [-1,-1,-1];
                [-1,-1,0];
                [-1,-1,1];

                [-1,0,-1];
                [-1,0,0];
                [-1,0,1];

                [-1,1,-1];
                [-1,1,0];
                [-1,1,1];
                ];


                Rank_voxel=nan(size(Nearborcor,1),size(datarotate_timecorrect_gf,4));

                        for nbi=1:size(Nearborcor,1)

                            %correct typo here 20231011
                            if ((sum([1:size(datarotate_timecorrect_gf,1)]==Nearborcor(nbi,1))>0)&&((sum([1:size(datarotate_timecorrect_gf,2)]==Nearborcor(nbi,2))>0))&&(sum([1:size(datarotate_timecorrect_gf,3)]==Nearborcor(nbi,3))>0))

                                comptrace=reshape(Rank(Nearborcor(nbi,1),Nearborcor(nbi,2),Nearborcor(nbi,3),:),[],1);
                                
                                Rank_voxel(nbi,:)=comptrace;
                            end

                        end

                        KCC_local(xi,yi,zi)=(12*sum(mean(Rank_voxel,'omitnan').^2)/(size(Rank,4)^3-size(Rank,4)))-3*((size(Rank,4)+1)/(size(Rank,4)-1));
                        
           end
        end
    end
    save(strcat('\\fenglab03\yiyun\20240820 fUS\n',S(si).name(1:strfind(S(si).name,'re')-1),'kcc.mat'),'Rank','KCC_local');

end

%%
