function PickedBlocks = selectBlocks(SlopeType,SlopeRate)
    % e.g. PickedBlocks = selectBlocks('cont',1);
    addpath(genpath('Data'))
    load('Data/IDs.mat')
    
    fprintf('Identify block trials... ')
    tic;
    for animal = 1:length(IDs)
        if exist(['Data_' IDs{animal} '.mat'])
            load(['Data_' IDs{animal} '.mat'])
        else
            continue
        end
        for session = 1:length(Data)
            if ~isfield(Data{session},'Trials')
                continue
            else
            end
            if ~isfield(Data{session}.Trials,'Start')
                continue
            else
            end
            if length(Data{session}.Trials.Start)<=16
                continue
            else
            end
            FindStart = 1;
            FindEnd = 0;
            countHR = 0;
            countLR = 0;
            countCut = 0;
            BlockStarts{animal}{session} = [];
            BlockEnds{animal}{session} = [];
            for trial = 1:length(Data{session}.Trials.Start)
                countCut = countCut + 1;
                if FindStart==1 && countCut<100
                    if Data{session}.Trials.HRreward(trial,1)~=0
                        countHR = countHR + 1;
                    else
                        countHR = 0;
                    end
                    if countHR>=4
                        BlockStarts{animal}{session} = [BlockStarts{animal}{session}; trial-3];
                        FindStart = 0;
                        FindEnd = 1;
                        countHR = 0;
                    else
                    end
                    countLR = 0;
                elseif FindEnd==1 || countCut>=100
                    if countCut>=100
                        BlockEnds{animal}{session} = [BlockEnds{animal}{session}; trial];
                        FindStart = 1;
                        FindEnd = 0;
                        countLR = 0;
                        countCut = 0;
                    else
                        if Data{session}.Trials.LRreward(trial,1)==1
                            countLR = countLR + 1;
                        else
                            countLR = 0;
                        end
                        if countLR>=6
                            nextTrial = 1;
                            if trial+nextTrial>length(Data{session}.Trials.LRreward)
                            else
                                while Data{session}.Trials.LRreward(trial+nextTrial,1)==1
                                    nextTrial = nextTrial + 1;
                                    if trial+nextTrial>length(Data{session}.Trials.LRreward)
                                        break
                                    else
                                    end
                                end
                            end
                            BlockEnds{animal}{session} = [BlockEnds{animal}{session}; trial+nextTrial-1];
                            FindStart = 1;
                            FindEnd = 0;
                            countLR = 0;
                            countCut = 0;
                        else
                        end
                    end
                    countHR = 0;
                else
                    error('What to find???')
                end
            end
            if isempty(BlockEnds{animal}{session}) || isempty(BlockStarts{animal}{session})
                continue
            else
            end
            checkcount = 0;
            while length(BlockEnds{animal}{session})>length(BlockStarts{animal}{session}) && checkcount<=3
                for block = 1:length(BlockStarts{animal}{session})
                    if BlockEnds{animal}{session}(block)-BlockStarts{animal}{session}(block)<0
                        BlockEnds{animal}{session}(block) = NaN;
                        break
                    else
                    end
                end
                BlockEnds{animal}{session}(isnan(BlockEnds{animal}{session})) = [];
                checkcount = checkcount + 1;
            end
            if checkcount>3
                continue
            else
            end
            if length(BlockEnds{animal}{session})==length(BlockStarts{animal}{session})
                for block = 1:length(BlockStarts{animal}{session})
                    if BlockEnds{animal}{session}(block)-BlockStarts{animal}{session}(block)<0
                        BlockEnds{animal}{session}(block) = NaN;
                        break
                    else
                    end
                end
                BlockEnds{animal}{session}(isnan(BlockEnds{animal}{session})) = [];
            else
            end
            Data{session}.Trials.BlockStart = BlockStarts{animal}{session};
            Data{session}.Trials.BlockEnd = BlockEnds{animal}{session};
        end
    end
    toc;

    fprintf('Identify block types... ')
    tic;
    for animal = 1:length(IDs)
        if exist(['Data_' IDs{animal} '.mat'])
            load(['Data_' IDs{animal} '.mat'])
        else
            continue
        end
        for session = 1:length(Data)
            if ~isfield(Data{session},'Trials')
                continue
            else
            end
            if ~isfield(Data{session}.Trials,'Start')
                continue
            else
            end
            if length(Data{session}.Trials.Start)<=20
                continue
            else
            end
            if length(Data{session}.Trials.BlockEnd)>length(Data{session}.Trials.BlockStart) || isempty(Data{session}.Trials.BlockEnd)
                continue
            else
            end
            for block = 1:length(Data{session}.Trials.BlockEnd)
                if Data{session}.Trials.BlockEnd(block)<=Data{session}.Trials.BlockStart(block)
                    continue
                else
                end
                Block{animal}{session}{block} = identifyBlock(animal,session,Data{session}.Trials.BlockStart(block),Data{session}.Trials.BlockEnd(block));
            end
        end
    end
    toc;
    
    fprintf('Selecting blocks... ')
    tic;
    c = 0;
    for animal = 1:length(IDs)
        if exist(['Data_' IDs{animal} '.mat'])
            load(['Data_' IDs{animal} '.mat'])
        else
            continue
        end
        for session = 1:length(Data)
            if ~isfield(Data{session},'Trials')
                continue
            else
            end
            if ~isfield(Data{session}.Trials,'Start')
                continue
            else
            end
            if length(Data{session}.Trials.Start)<=20
                continue
            else
            end
            if length(Data{session}.Trials.BlockEnd)>length(Data{session}.Trials.BlockStart) || isempty(Data{session}.Trials.BlockEnd)
                continue
            else
            end
            for block = 1:length(Block{animal}{session})
                if Data{session}.Trials.BlockEnd(block)<=Data{session}.Trials.BlockStart(block)
                    continue
                else
                end
                if contains(Block{animal}{session}{block}.SlopeType,SlopeType) && Block{animal}{session}{block}.SlopeRate==SlopeRate
                    c = c + 1;
                    PickedBlocks.Features(1,c) = Block{animal}{session}{block}.HRsigmoid_coefficient(2);
                    PickedBlocks.Features(2,c) = Block{animal}{session}{block}.Onset;
                    PickedBlocks.Features(3,c) = Block{animal}{session}{block}.MeanValue;
                    PickedBlocks.Features(4,c) = Block{animal}{session}{block}.MidPoit;
                    PickedBlocks.Features(5,c) = Block{animal}{session}{block}.Offnset;
                    PickedBlocks.Features(6,c) = Block{animal}{session}{block}.MaxHRpress;
                    PickedBlocks.Features(7,c) = length(Block{animal}{session}{block}.Trials); % Block length
                    PickedBlocks.Features(8,c) = nanmedian(Block{animal}{session}{block}.Initiation(2:end)-Block{animal}{session}{block}.Rewarded(1:end-1));% Median inter-trial interval
                    PickedBlocks.Features(9,c) = Block{animal}{session}{block}.MissRate + Block{animal}{session}{block}.AbortRate;% Miss+Abort rate
                    PickedBlocks.AnimalSessionBlock(c,:) = [animal session block];
                else
                end
            end
        end
    end
    PickedBlocks.Picked = nanmax(abs((PickedBlocks.Features-nanmean(PickedBlocks.Features,2))./nanstd(PickedBlocks.Features,[],2)));
    PickedBlocks.Picked(PickedBlocks.Picked<1.0365) = 1; % 70% of data
    PickedBlocks.Picked(PickedBlocks.Picked>=1.0365) = 0;
    Traces = [];
    for blk = 1:length(PickedBlocks.Picked)
        if PickedBlocks.Picked(blk)==1
            Traces = [Traces; Block{PickedBlocks.AnimalSessionBlock(blk,1)}{PickedBlocks.AnimalSessionBlock(blk,2)}{PickedBlocks.AnimalSessionBlock(blk,3)}.HRsigmoid];
        else
        end
    end
    toc;
end

function Block = identifyBlock(Animal,Session,StartTrial,EndTrial)
    load('Data/IDs.mat')
    load(['Data_' IDs{Animal} '.mat'])
    
    trials = StartTrial:EndTrial;
    Block.File = Data{Session}.File;
    Block.Trials = trials';
    Block.Start = Data{Session}.Trials.Start(trials);
    Block.End = Data{Session}.Trials.End(trials);
    Block.HRchoice = Data{Session}.Trials.HRchoice(trials);
    Block.LRchoice = Data{Session}.Trials.LRchoice(trials);
    Block.HRpress = Data{Session}.Trials.HRpress(trials);
    Block.LRpress = Data{Session}.Trials.LRpress(trials);
    Block.HRreward = Data{Session}.Trials.HRreward(trials);
    Block.LRreward = Data{Session}.Trials.LRreward(trials);
    Block.OptoON = Data{Session}.Trials.OptoON(trials);
    Block.OptoOFF = Data{Session}.Trials.OptoOFF(trials);
    Block.Initiation = Data{Session}.Trials.Initiation(trials);
    Block.Rewarded = Data{Session}.Trials.Rewarded(trials);
    Block.Slope = Data{Session}.Trials.Slope(trials);
    opto = Block.OptoON;
    opto(isnan(opto)) = [];
    if length(opto)>5
        Block.OptoType = Data{Session}.Info.OptoType;
        if contains(Data{Session}.Info.Note,'laser') || contains(Data{Session}.Info.Note,'Laser') || contains(Data{Session}.Info.Note,'drop') || contains(Data{Session}.Info.Note,'fiber')  || contains(Data{Session}.Info.Note,'Fiber')
            Block.OptoType = 'NonOpto';
        else
        end
    else
        if contains(Data{Session}.Info.OptoType,'22Q11')
            Block.OptoType = '22Q11';
        else
            Block.OptoType = 'NonOpto';
        end
    end
    design = Data{Session}.Info.OptoDesign;
    if contains(design,'PR')
        des1 = 'PR';
    elseif contains(design,'FR')
        des1 = 'FR';
    else
        des1 = 'All';
    end
    if contains(design,'post-')
        des2 = 'post-reward';
    elseif contains(design,'pre-')
        des2 = 'pre-reward';
    else
        des2 = 'NA';
    end
    if contains(Block.OptoType,'NonOpto') || contains(Block.OptoType,'22Q11')
        Block.OptoDesign = 'NA';
    else
        Block.OptoDesign = [des1 ' ' des2];
    end
    if length(opto)>5
        if contains(Data{Session}.Info.Note,'laser') || contains(Data{Session}.Info.Note,'Laser') || contains(Data{Session}.Info.Note,'drop') || contains(Data{Session}.Info.Note,'fiber')  || contains(Data{Session}.Info.Note,'Fiber')
            Block.OptoType = 'NA';
        else
        end
    else
    end
    
    [Block.HRsigmoid_coefficient,model] = fitSigmoid5(Block.HRreward);
    Block.HRsigmoid = model(Block.HRsigmoid_coefficient,1:100);
    finescaleModel = model(Block.HRsigmoid_coefficient,0.001:0.001:100);
    [~,Block.Onset] = min(abs(finescaleModel-0.9));
    Block.Onset = Block.Onset/1000;
    [~,Block.IndiffPoint] = min(abs(finescaleModel-0.5));
    Block.IndiffPoint = Block.IndiffPoint/1000;
    [~,Block.Offnset] = min(abs(finescaleModel-0.1));
    Block.Offnset = Block.Offnset/1000;
    Block.MeanValue = (nansum(Block.HRreward)+nansum(Block.LRreward))./(nansum(Block.HRpress)+nansum(Block.LRpress));
    Block.MaxHRpress = nanmax(Block.HRpress);
    for trial = 1:length(Block.Trials)
        if Block.HRreward(trial,1)+Block.LRreward(trial,1)==0
            if Block.HRchoice(trial,1)+Block.LRchoice(trial,1)==2
                Block.Abort(trial,1) = 1;
            else
                Block.Miss(trial,1) = 1;
            end
        else
            Block.Miss(trial,1) = 0;
            Block.Abort(trial,1) = 0;
        end
    end
    Block.MissRate = sum(Block.Miss)./length(Block.Miss);
    Block.AbortRate = sum(Block.Abort)./length(Block.Abort);
    
    Block.HRpattern = Block.HRpress(Block.HRreward~=0);
    pattern = Block.HRpattern(4:end);
    delta = unique(pattern(2:end)-pattern(1:end-1));
    if length(delta)==1 && delta==1
        Block.SlopeType = 'cont';
        Block.SlopeRate = 1;
    elseif length(delta)==1 && delta==2
        Block.SlopeType = 'cont';
        Block.SlopeRate = 2;
    elseif length(delta)==1 && delta==3
        Block.SlopeType = 'cont';
        Block.SlopeRate = 3;
    elseif length(delta)==1 && delta==0
        Block.SlopeType = 'cont';
        Block.SlopeRate = 0;
    elseif length(delta)==2 && max(delta)==1
        Block.SlopeType = 'disc';
        Block.SlopeRate = 0.5;
    elseif length(delta)==2 && max(delta)==2
        Block.SlopeType = 'disc';
        Block.SlopeRate = 1;
    elseif length(delta)==2 && max(delta)==4
        if ~isempty(strfind([pattern(2:end)-pattern(1:end-1)]',[0 0 0 4]))
            Block.SlopeType = '4td';
            Block.SlopeRate = 1;
        elseif ~isempty(strfind([pattern(2:end)-pattern(1:end-1)]',[0 4]))
            Block.SlopeType = 'disc';
            Block.SlopeRate = 2;
        else
            Block.SlopeType = 'unclassified';
            Block.SlopeRate = NaN;
        end
    elseif length(delta)==2 && max(delta)==6
        Block.SlopeType = 'disc';
        Block.SlopeRate = 3;
    else
        Block.SlopeType = 'unclassified';
        Block.SlopeRate = NaN;
    end
    
    if ~isnan(Block.SlopeRate)
        if contains(Block.SlopeType,'cont') || contains(Block.SlopeType,'disc')
            [Block.HRpolicy,Block.HRpolicy_coefficient,~] = estimatePolicyV2(IDs{Animal},Session,StartTrial,EndTrial,Block.SlopeRate,Block.SlopeType);
        else
            Block.HRpolicy = NaN;
            Block.HRpolicy_coefficient = NaN;
        end
    else
        Block.HRpolicy = NaN;
        Block.HRpolicy_coefficient = NaN;
    end
end

function [k,model] = fitSigmoid5(BinaryData)
    type = 'Normal';
    model = @(k,x)  (1-cdf(type,x,k(1),k(2)));
    k = nlinfit([1:length(BinaryData)]',smoothdata(BinaryData,'gaussian',6,'omitnan'),model,[length(BinaryData)/2 1],'Options',statset('FunValCheck','off','MaxIter',20000,'WgtFun','logistic'));
    k(2) = k(2)/3;
    Y_model = model(k,[1:length(BinaryData)]');
end