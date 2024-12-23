function Info = parseLeverPress(fn,animalID)
%Input:
%       fn - txt file with or without path. Required
%       ID - animal ID. Not required 
%Output: Info struct with the following fields:
%       ID - animal ID. If not in input the ID is left blank
%       LRSide - side of low reward lever (left/right)
%       LRrequest - num of presses required for low reward
%       HRSide - side of high reward (left/right)
%       HRtrialsToIncrease - number of consecutive trials required to start
%           increasing number of presses required for high reward
%       HRrequest - number of trials for high reward. Broken into cell array
%           for each increasing block
%       VolumeRatio - ratio of high reward amount to low reward amount
%       TrialList - cell array of trial order. Listed as LR, HR, and Miss
%           (unrewarded)
%       Trials.LRchoice - binary array of trials in which the animal
%           choce the LR side first
%       Trials.HRchoice - binary array of trials in which the animal chose
%           HR side first
%       Trials.LRpresses -  array of number of presses on the LR side for
%           each trial
%       Trials.HRpresses  -  array of number of presses on the HR side for
%           each trial
%       Trials.LRreward - binary array of trials in which the animal
%           recieved low reward
%       Trials.HRreward  - binary array of trials in which the animal
%           recieved high reward
%       Trials.Opto- binary array of trials with laser on
%       Trials.BlockNumList - array of which block each trial belongs to
%       Trials.BlockStart - binary array of the first trial in each block
%           (includes first trial of session)
%       Trials.BlockEnd - binary array of final trial in each block
%           (includes final trial of session)

[behavior] = textread(fn,'%s',-1,'delimiter','\t');

if ~exist('animalID')
    animalID = ' ';
end

get_text = @(x,y) (cellfun(@(c) contains(c,y),x));

trial_starts = get_text(behavior,'trial available');
%cellfun(@(c) contains(c,'trial available'),behavior);
trial_ind = find(trial_starts==1);
trial_ind = [trial_ind; size(trial_starts,1)];
num_trials = size(trial_ind,1)-1;

FRside = [];PRside = []; RewardAmt = {}; TrialList = {}; tt = 1;
for t = 1:num_trials
    trial_text = behavior(trial_ind(t):trial_ind(t+1)-1);
    
    removeChk = get_text(trial_text,'removed');
    if sum(removeChk) > 0
        continue;
    end
    
    initChk = get_text(trial_text,'initiation');
    if sum(initChk) == 0
        continue;
    end
    
    blockChk = get_text(trial_text,'blocknumber');
    BlockTxt = trial_text(blockChk);
    BlockNum(tt) = str2double(BlockTxt{1}(end));
    
    completeChk = get_text(trial_text,'reward delivered');
    if sum(completeChk) > 0
        RewardType = trial_text(get_text(trial_text,'reward available'));
        RewardAmtIdx = trial_text(get_text(trial_text,'amt'));
        if contains(RewardType,'fixed')
            TrialType = 'LR';
            RewardAmt{tt} = RewardAmtIdx{1}(end);
        elseif contains(RewardType,'progressive')
            TrialType = 'HR';
            RewardAmt{tt} = RewardAmtIdx{1}(end);
        end
    else
        TrialType = 'Miss';
        RewardAmt{tt} = [];
    end
          
    TrialList{tt} = TrialType;
    
    pokedChk = get_text(trial_text,'poked');
    if sum(pokedChk) == 0
        numPresses(tt) = 0;
        numPresses_left(tt) = 0;
        numPresses_right(tt) = 0;
        
        leftchoice(tt) = 0;
        rightchoice(tt) = 0;
        Side = 'NaN';

    else
        pokedtext = trial_text(pokedChk);
        numPresses(tt) = numel(pokedtext);
        choicepoke = pokedtext{1};
        pokedChk_left = get_text(pokedtext,'left');
        pokedChk_right = get_text(pokedtext,'right');
        
        numPresses_left(tt) = sum(pokedChk_left);
        numPresses_right(tt) = sum(pokedChk_right);
                
        if contains(choicepoke,'left')
            leftchoice(tt) = 1;
            rightchoice(tt) = 0;
        elseif contains(choicepoke,'right')
            leftchoice(tt) = 0;
            rightchoice(tt) = 1;
        end
        
        if (numPresses_left(tt) == numPresses(tt))
            Side = 'Left';
        elseif (numPresses_right(tt) == numPresses(tt))
            Side = 'Right';
        else
            Side = 'NaN';
        end
        
    end
        
    
    switch TrialType
        case 'LR'
            if isempty(FRside)
                FRside = Side;
            end
        case 'HR'
            if isempty(PRside)
                PRside = Side;
            end
        case 'Miss'
    end
    
    LaserOnChk = get_text(trial_text,'laser on');
    if sum(LaserOnChk) == 0
        LaserTrial(tt) = 0;
    elseif sum(LaserOnChk) > 0
        LaserTrial(tt) = 1;
    end
    
    tt = tt+1;
    
end
TrialList = TrialList';
numPresses = numPresses';
RewardAmt = str2double(RewardAmt)';
BlockNum = BlockNum'+1;
BlockDiff = diff(BlockNum);
BlockSwitch = find(BlockDiff);


if strcmp(FRside,'Right')
    LRchoice = rightchoice';
    HRchoice = leftchoice';
    LRpresses = numPresses_right';
    HRpresses = numPresses_left';
else
    LRchoice = leftchoice';
    HRchoice = rightchoice';
    LRpresses = numPresses_left';
    HRpresses = numPresses_right';
end
    
LRreward = (strcmp(TrialList,'L:R') & ~isnan(RewardAmt));
HRreward = (strcmp(TrialList,'HR') & ~isnan(RewardAmt));

PRidx = find(get_text(TrialList,'HR'));
FRidx = find(get_text(TrialList,'LR'));
FRrequest = numPresses(FRidx(1));
VolumeRatio = RewardAmt(PRidx(1))/RewardAmt(FRidx(1));

tempTrialList = [TrialList; 'Empty'];
for f = 1:numel(BlockSwitch)
    reset = BlockSwitch(f);
    for ff = reset:numel(tempTrialList)
        testtrial1 = tempTrialList{ff};
        testtrial2 = tempTrialList{ff+1};
        
        if strcmp(testtrial2,testtrial1)
            continue;
        else
            BlockEndIdx(f) = ff;
            break;
        end
    end
end

PRidxx = [0;PRidx];
PRpresses = numPresses(PRidx);
prcounter = 0;
InBlock = 0;
pp = 1;
for p = 1:numel(PRidx)
    testPress = PRpresses(p);
    if (testPress > 1 && InBlock == 0)
        lastbeforeincrease = p;
        for x = lastbeforeincrease:-1:2
            testtrial1 = PRidxx(x);
            testtrial2 = PRidxx(x-1);
            
            prcounter = prcounter+1;
            if (testtrial1-testtrial2) == 1
            else
                break;
            end
        end
        InBlock = 1;
        PRTrialsToIncreaseByBlock(pp) = prcounter;
        PRstartIncrease(pp) = lastbeforeincrease;
        prcounter = 0;
        pp = pp+1;
    end
    
    if testPress == 1
        InBlock = 0;
    end
end
PRTrialsToIncrease = min(PRTrialsToIncreaseByBlock);
PRlow = min(PRpresses);

PRpresstemp = [PRpresses;PRlow];
PRlowIdx = find(PRlow == PRpresstemp);
for y = 1:numel(PRstartIncrease)
    tempIncrease = PRstartIncrease(y);
    tempStart = tempIncrease - PRTrialsToIncreaseByBlock(y);
    tempEnd = find(PRlowIdx > tempIncrease,1,'first');
    
    BlockStartIdx(y) = PRidx(tempStart);
    
    PRrequest{y} = PRpresstemp(tempStart:PRlowIdx(tempEnd)-1);
end

BlockVals = zeros(numel(TrialList),1);
BlockStart = BlockVals;
BlockStart(BlockStartIdx) = 1;
BlockEnd = BlockVals;
BlockEnd(BlockEndIdx) = 1;

BlockList = zeros(numel(TrialList),1);
blockcounter = 1;
for bc = 1:numel(TrialList)
    if (blockcounter > numel(BlockStartIdx)) || (blockcounter > numel(BlockEndIdx))
        break;
    end
    
    tempBlockStart = BlockStartIdx(blockcounter);
    tempBlockEnd = BlockEndIdx(blockcounter);
    
    if (bc >= tempBlockStart) && (bc <=tempBlockEnd)
        BlockList(bc) = blockcounter;
    end
    if bc == tempBlockEnd
        blockcounter = blockcounter+1;
    end
end

Info.ID = animalID;
Info.LRSide = FRside;
Info.LRrequest = FRrequest;
Info.HRSide = PRside;
Info.HRtrialsToIncrease = PRTrialsToIncrease;
Info.HRrequest = PRrequest;
Info.VolumeRatio = VolumeRatio;
Info.TrialList = TrialList;
Info.Trials.LRchoice = LRchoice;
Info.Trials.HRchoice = HRchoice;
Info.Trials.LRpresses = LRpresses;
Info.Trials.HRpresses = HRpresses;
Info.Trials.LRreward = LRreward;
Info.Trials.HRreward = HRreward;
Info.Trials.Opto = LaserTrial';
Info.Trials.BlockNumList = BlockList;
Info.Trials.BlockStart = BlockStart;
Info.Trials.BlockEnd = BlockEnd;
end