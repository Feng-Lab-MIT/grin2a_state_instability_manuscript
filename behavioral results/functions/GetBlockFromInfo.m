
%this function if for getting block information from info got from
%getInfoFromTxt_0712, the difference with the previous version is that I
%included all trials in each block, which means we can also see
%inter-block-interval trials.  

% The reason why I used this new block identification method is that the
% previous one will not include the step function blocks. 
% with this new code I can 1. plot new blocks and pick step function block
% manually 2. Get interblock interval, which will be important for studying
% the cognitive flexibility/belief update

function data = GetBlockFromInfo(info)


if ~isempty(info)&~isempty(info.BlockStart)
    if length(info.BlockStart)>length(info.BlockEnd)
      
    info.BlockEnd(end+1)=length(info.LRchoice);
    end
    for i=1:length(info.BlockStart)
        
    block{i}=info.BlockStart(i):info.BlockEnd(i);
    end



%get block  info from session info

clear i
for i=1:length(block)
    data{i}.ID=info.ID;
    data{i}.File=info.File;
    data{i}.RewardedTimes=info.RewardedTimes(block{i});
    data{i}.InitiationTimes = info.InitiationTimes(block{i})-nanmin(info.InitiationTimes(block{i}));
    data{i}.LRrequest = info.LRrequest;
%     data{i}.HRpattern =info.HRpattern{i};
    data{i}.VolumeRatio = info.VolumeRatio;
    data{i}.LRchoice = info.LRchoice(block{i});
    data{i}.HRchoice = info.HRchoice(block{i});
    data{i}.LRpress = info.LRpress(block{i});
    data{i}.HRpress = info.HRpress(block{i});
    data{i}.LRreward = info.LRreward(block{i});
    data{i}.HRreward = info.HRreward(block{i});
    data{i}.LaserON = info.LaserON(block{i});
    data{i}.Manipulation = info.Manipulation;
    %data{i}.SlopeType = info.SlopeType{i};
   data{i}.HRrequest = info.HRpress(block{i});
    data{i}.PressGoals = info.PressGoals(block{i});
    data{i}.LaserON = info.LaserON(block{i});
    

    for j=2:length(data{i}.HRrequest)
        if data{i}.HRreward(j)>0
            data{i}.HRrequest(j)=data{i}.HRpress(j);
        else data{i}.HRrequest(j)=data{i}.HRrequest(j-1);
        end
    end
    % find the acurate blockEnd
    t=findConSame(data{i}.LRreward,1,6);
    if t>4
        blockend=min(t)+5;
            else
        blockend=max(find(data{i}.HRrequest==max(data{i}.HRrequest)));
    end
        data{i}.HRrequest=data{i}.HRrequest(1:blockend);
        data{i}.HRchoice=data{i}.HRchoice(1:blockend);
        data{i}.LRchoice = data{i}.LRchoice(1:blockend);
        data{i}.LRpress = data{i}.LRpress(1:blockend);
        data{i}.LRreward = data{i}.LRreward(1:blockend);
        data{i}.HRreward = data{i}.HRreward(1:blockend);
        data{i}.HRpress = data{i}.HRpress(1:blockend);
        data{i}.LaserON = data{i}.LaserON(1:blockend);
        data{i}.PressGoals = data{i}.PressGoals(1:blockend);

        
            
 
     
    
%     temp=max(find(data{i}.HRrequest==1));
%     data{i}.HRrequest=data{i}.HRrequest(temp:end);
%     data{i}.HRchoice=data{i}.HRchoice(temp:end);
%     clear temp
    if ~isempty(data{i}.HRrequest)
    [data{i}.Policy.X_HRrequest,data{i}.Policy.Y_pHR] = estimatePolicy(data{i}.HRrequest,data{i}.HRchoice);
    else data{i}.Policy.X_HRrequest=[];data{i}.Policy.Y_pHR=[];
    end
    
end
else
    data=[];
end


end

