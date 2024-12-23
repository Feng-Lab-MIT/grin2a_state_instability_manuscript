%1. get raw trace from each slices
%2. get selected time window <<<
%3. interpolate => filter  <<<
%4. rotation and align


% filtering << and baseline remove


%to do: select time window <<< 
%

%=============================

%datascan = openScan('Z:\fUSdata\GRIN2A\scans\sub-9854_1\ses-Session_2023-2-25\4Dscan_1_9854_1_fus3D.source.scan');

S=dir(fullfile('Z:\fUSdata\fus_rawdata_for_connect_no_dex\raw_data\','*.h5'));

Seltimewin=repmat([800 1400],numel(S),1);
% Seltimewin(1,:)=[0 600];
% Seltimewin([7,12,13,14],:)=repmat([600 1200],4,1);


%%
% for checking where to cut %20-100
% figure()
% for si=1:numel(S)
% %
% rawfilename=strcat('Z:\fUSdata\fus_rawdata_for_connect_no_dex\raw_data\',S(si).name);
% 
% dataori = h5read(rawfilename,'/Acquisition');
% data2 = permute(reshape(dataori, fliplr(size(dataori))), [3 1 2 4]);
% 
% 
% datapreroate=zeros(size(data2,2),size(data2,3),size(data2,1),size(data2,4));
% for ti=1:size(data2,4)
%     datapreroate(:,:,:,ti) = shiftdim(data2(:,:,:,ti),1);
% end
% 
% 
% for sli=1:7
%     figure(sli)
%     subplot(6,6,si)
%     imagesc(reshape(datapreroate(:,sli,:,1),size(datapreroate,1),[])');
% end
% 
% 
% end


%%

for si=1:numel(S) %7,11,

%% Import data    
%rawfilename='Z:\fUSdata\fus_rawdata_for_connect\export_data_98547.h5';
%rawfilename=strcat('Z:\fUSdata\fus_rawdata_for_connect\',S(si).name);
rawfilename=strcat('Z:\fUSdata\fus_rawdata_for_connect_no_dex\raw_data\',S(si).name);

info = h5info(rawfilename);
timeold = h5read(rawfilename,'/TimeOld');
time = h5read(rawfilename,'/Time');

dataori = h5read(rawfilename,'/Acquisition');
data2 = permute(reshape(dataori, fliplr(size(dataori))), [3 1 2 4]);

braintolab = h5read(rawfilename,'/BrainToLab');

braintolab=reshape(braintolab,4,[]);
braintolabinv= inv(braintolab); %rotation matrix

voxeltoprobe = h5read(rawfilename,'/VoxelsToProbe');
probetolab = h5read(rawfilename,'/ProbeToLab');

%% Import time points
initime=Seltimewin(si,1);
fintime=Seltimewin(si,2);

%% select window
timeoldresort=sort(timeold);

leftwind=max([find(timeoldresort>=initime,1,'first'),find(timeoldresort>=timeold(2),1,'first')]);
rightwind=min([find(timeoldresort<=fintime,1,'last'),find(timeoldresort<=timeold(end-size(data2,3)+1),1,'last')]);

leftwindt=find(time>=timeoldresort(leftwind),1,'first');
rightwindt=find(time<=timeoldresort(rightwind),1,'last');

%% filter

Fs=1/median(diff(timeoldresort));
d = designfilt('bandpassiir', ...       % Response type
       'DesignMethod','butter', ...  % Design method
       'FilterOrder',4,...  %
       'HalfPowerFrequency1',0.002,...  %0.0008  %change this to higher
       'HalfPowerFrequency2',0.12,...  %change from 0.15 to 0.12 23/10/11
       'SampleRate',Fs);   


%% get coordinates

datapreroate=zeros(size(data2,2),size(data2,3),size(data2,1),size(data2,4));
for ti=1:size(data2,4)
    datapreroate(:,:,:,ti) = shiftdim(data2(:,:,:,ti),1);
end

datapreroate_timecorrect=zeros(size(datapreroate,1),size(datapreroate,2),size(datapreroate,3),rightwindt-leftwindt+1);

%startxi=20;%%%%% start from 60
%endxi=100;%%%%% start from 60

xyzcor=zeros(size(datapreroate,1)*size(datapreroate,2)*size(datapreroate,3),3);
olddata_cor=zeros(size(datapreroate,1)*size(datapreroate,2)*size(datapreroate,3),4);
olddata_timeser=zeros(size(datapreroate,1)*size(datapreroate,2)*size(datapreroate,3),rightwind-leftwind+1);
filtolddata_timeser=zeros(size(datapreroate,1)*size(datapreroate,2)*size(datapreroate,3),rightwindt-leftwindt+1);
transformeddata_cor=zeros(size(datapreroate,1)*size(datapreroate,2)*size(datapreroate,3),4);

k=1;



for xi=1:size(datapreroate,1) %startxi:endxi%
    for yi=1:size(datapreroate,2)
        for zi=1:size(datapreroate,3)
            
            xyzcor(k,:)=[xi,yi,zi];

            olddata_cor(k,:)=[xi,yi,zi,1]*reshape(voxeltoprobe,4,[])*reshape(probetolab,4,[]);
            
            %=================================================
            %===interp to get new data points for time========
            % time:downsampled (one for whole volume)
            % timeold: timepoint of each frames
            %=================================================
            timeoldresort_cut=timeoldresort(leftwind:rightwind); %this only consider interested window
            olddata_timeser(k,:)=interp1(timeold(yi:size(datapreroate,2):length(timeold)),reshape(datapreroate(xi,yi,zi,:),[],1),timeoldresort_cut); 
            
            %=================================================
            %===filtering========
            %=================================================
            filtolddatatimeser=filtfilt(d,olddata_timeser(k,:));

            %=================================================
            %===remove extreme value========
            %=================================================
                
            windw=2;%time before and after extreme value to be removed with (second)
            extremefold=2;%when signal larger than how much fold of std  
            timepoint_to_remove=[];
            extreme_data_time=find(filtolddatatimeser>=(extremefold*std(filtolddatatimeser)+mean(filtolddatatimeser)));

            for exti=1:length(extreme_data_time)
                onset=find(timeoldresort_cut<(timeoldresort_cut(extreme_data_time(exti))-windw),1,'last'); 
                offset=find(timeoldresort_cut>(timeoldresort_cut(extreme_data_time(exti))+windw),1,'first');
            
                timepoint_to_remove=[timepoint_to_remove,onset:offset];
            end
            
            timepoint_to_remove=unique(timepoint_to_remove);
            %timepoint_to_remove(timepoint_to_remove<leftwind)=[];
            %timepoint_to_remove(timepoint_to_remove>min([rightwind,length(filtolddatatimeser)]))=[];
            timepoint_to_remove(timepoint_to_remove<1)=[];
            timepoint_to_remove(timepoint_to_remove>length(filtolddatatimeser))=[];
           
            
            timepoint_to_use=1:length(filtolddatatimeser);
            
            timepoint_to_use(timepoint_to_remove)=[];
            
            %fill removed time points with new data
            filtolddatatimeser2=interp1(timeoldresort_cut(timepoint_to_use),filtolddatatimeser(timepoint_to_use),timeoldresort_cut(1:length(filtolddatatimeser)),'pchip'); %this is wrong

            %=================================================
            %===downsample========
            %=================================================
            %downsample
            filtolddata_timeser(k,:)=interp1(timeoldresort_cut,filtolddatatimeser2,time(leftwindt:rightwindt));

            %=================================================
            %===rotate to align to common atlas========
            %=================================================            

            transformeddata_cor(k,:)=[xi,yi,zi,1]*reshape(voxeltoprobe,4,[])*reshape(probetolab,4,[])*braintolabinv; %rotation
            
            datapreroate_timecorrect(xi,yi,zi,:)=filtolddata_timeser(k,:);
            
            
            k=k+1;
        end
    end
end

%save(strcat('D:\20230515 fUS no dex\local_connectivity_outliner_remove\n',rawfilename(64:strfind(rawfilename,'.')-1),'temfiltnrotat1011.mat'),'olddata_cor','datapreroate_timecorrect','transformeddata_cor','filtolddata_timeser');

%% resamaple first (rotate and resample)

%inteprolate

datarotate_timecorrect=zeros(length([-1.2:0.03:1.2]),length([-0.6:0.03:0.6]),length([-0.8:0.03:1.4]),size(datapreroate_timecorrect,4));

for timei=1:size(datapreroate_timecorrect,4)
    Frot = scatteredInterpolant(transformeddata_cor(:,1),transformeddata_cor(:,2),transformeddata_cor(:,3),filtolddata_timeser(:,timei));
    
    [xq,yq,zq] = meshgrid(-0.6:0.03:0.6,-1.2:0.03:1.2,-0.8:0.03:1.4);  % you may need to adjust here
    datarotate_timecorrect(:,:,:,timei) = Frot(xq,yq,zq);

    clearvars Frot 
end

save(strcat('D:\20230515 fUS no dex\local_connectivity_outliner_remove\n',rawfilename(64:strfind(rawfilename,'.')-1),'temfiltnrotat1011.mat'),'olddata_cor','datapreroate_timecorrect','transformeddata_cor','filtolddata_timeser','datarotate_timecorrect');
%save(strcat('D:\20230515 fUS no dex\local_connectivity_outliner_remove\n',rawfilename(64:strfind(rawfilename,'.')-1),'globalcon0717.mat'),'Globalconnect','Fglob','vq_glob','olddata_cor','datapreroate_timecorrect','transformeddata_cor','filtolddata_timeser');

clearvars -except S Seltimewin datarotate_timecorrect xyzcor rawfilename


%% spatial filter: gaussina filter...  
% 
%data after temporal filter: datapreroate_timecorrect

datarotate_timecorrect_gf=zeros(size(datarotate_timecorrect,1),size(datarotate_timecorrect,2),size(datarotate_timecorrect,3),size(datarotate_timecorrect,4));

%need to resample to get homogeneous sample before 3d filtering... 

for timei=1:size(datarotate_timecorrect,4)


    datarotate_timecorrect_gf(:,:,:,timei)=imgaussfilt3(datarotate_timecorrect(:,:,:,timei),1);

%     figure()
%     for i=1:2:19
%     subplot(2,10,(i-1)/2+1)
%     imagesc(reshape(datarotate_timecorrect(:,i,:,timei),81,[]));
%     subplot(2,10,10+(i-1)/2+1)
%     imagesc(reshape(datarotate_timecorrect_gf(:,i,:,timei),81,[]));
%     end

end

save(strcat('D:\20230515 fUS no dex\local_connectivity_outliner_remove\n',rawfilename(64:strfind(rawfilename,'.')-1),'resamplengf1011.mat'),'datarotate_timecorrect_gf');
%save(strcat('D:\20230515 fUS no dex\local_connectivity_outliner_remove\n',rawfilename(64:strfind(rawfilename,'.')-1),'globalcon0717.mat'),'Globalconnect','Fglob','vq_glob','olddata_cor','datapreroate_timecorrect','transformeddata_cor','filtolddata_timeser');

clearvars -except S Seltimewin datarotate_timecorrect_gf xyzcor rawfilename


%% local connectivity

vq_local=zeros(size(datarotate_timecorrect_gf,1),size(datarotate_timecorrect_gf,2),size(datarotate_timecorrect_gf,3));


for xi=1:size(datarotate_timecorrect_gf,1)
    for yi=1:size(datarotate_timecorrect_gf,2)
        for zi=1:size(datarotate_timecorrect_gf,3)
    
    rawtrace=datarotate_timecorrect_gf(xi,yi,zi,:);
    

    Nearborcor=[xi,yi,zi]+[
%     [0,-1,-1];
%     [0,-1,0];
%     [0,-1,1];

    [0,0,1];
    [0,0,-1];

%     [0,1,-1];
%     [0,1,0];
%     [0,1,1];

%     [1,-1,-1];
%     [1,-1,0];
%     [1,-1,1];

    [1,0,-1];
    [1,0,0];
    [1,0,1];

%     [1,1,-1];
%     [1,1,0];
%     [1,1,1];

%     [-1,-1,-1];
%     [-1,-1,0];
%     [-1,-1,1];

    [-1,0,-1];
    [-1,0,0];
    [-1,0,1];

%     [-1,1,-1];
%     [-1,1,0];
%     [-1,1,1];
    ];

    
    Nearborxcor=nan(size(Nearborcor,1),1);

            for nbi=1:size(Nearborcor,1)
        
                %correct typo here 20231011
                if ((sum([1:size(datarotate_timecorrect_gf,1)]==Nearborcor(nbi,1))>0)&&((sum([1:size(datarotate_timecorrect_gf,2)]==Nearborcor(nbi,2))>0))&&(sum([1:size(datarotate_timecorrect_gf,3)]==Nearborcor(nbi,3))>0))
                  
                    comptrace=reshape(datarotate_timecorrect_gf(Nearborcor(nbi,1),Nearborcor(nbi,2),Nearborcor(nbi,3),:),[],1);
                    Coef=corrcoef(rawtrace(1:end-1),comptrace(1:end-1));
                    Nearborxcor(nbi)=Coef(2,1);
                end
        
            end
            
            vq_local(xi,yi,zi)=mean(Nearborxcor,'omitnan');

        end
    end
end


%%
save(strcat('D:\20230515 fUS no dex\local_connectivity_outliner_remove\n',rawfilename(64:strfind(rawfilename,'.')-1),'localcon_gf1011.mat'),'vq_local');

clearvars -except S Seltimewin

end