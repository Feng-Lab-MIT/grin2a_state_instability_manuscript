function Info = getInfoFromTxt(fn)

[behavior] = textread(fn,'%s',-1,'delimiter','\t');
[~,filename] = fileparts(fn);

filename_parts = split(filename,'_');
animalNameChk = cellfun(@(x) contains(x,'TT'),filename_parts);
%-------------------------------------------------------------------
animalName = filename_parts{animalNameChk};
%-------------------------------------------------------------------
get_text = @(x,y) (cellfun(@(c) contains(c,y),x));

trial_starts = get_text(behavior,'trial available');
%cellfun(@(c) contains(c,'trial available'),behavior);
trial_ind = find(trial_starts==1);
trial_ind = [trial_ind; size(trial_starts,1)];
num_trials = size(trial_ind,1)-1;

FRside = [];PRside = []; RewardAmt = {}; TrialList = {};
for t = 1:num_trials
    trial_text = behavior(trial_ind(t):trial_ind(t+1)-1);
    
    removeChk = get_text(trial_text,'removed');
    if sum(removeChk) > 0
        continue;
    end
    
    completeChk = get_text(trial_text,'reward delivered');
    if sum(completeChk) > 0
        RewardType = trial_text(get_text(trial_text,'reward available'));
        RewardAmtIdx = trial_text(get_text(trial_text,'Reward amt'));
        if contains(RewardType,'fixed')
            TrialType = 'FR';
            RewardAmt{t} = RewardAmtIdx{1}(end);
        elseif contains(RewardType,'progressive')
            TrialType = 'PR';
            RewardAmt{t} = RewardAmtIdx{1}(end);
        end
    else
        TrialType = 'Miss';
        RewardAmt{t} = [];
    end
          
    
%     if contains(RewardType,'fixed')
%         TrialType = 'FR';             
%          RewardAmt{t} =  RewardAmtIdx{1}(end);
%     elseif contains(RewardType,'progressive')
%         TrialType = 'PR';
%          RewardAmt{t} =  RewardAmtIdx{1}(end);
%     else
%         TrialType = 'Miss';
%         RewardAmt{t} = [];
%     end
    TrialList{t} = TrialType;
    
    if ~strcmp(TrialType,'Miss')
        pokedIdx = find(get_text(trial_text,'poked'));
        pokedtext = trial_text(pokedIdx(1));
        if contains(pokedtext,'right')
            Side = 'Right';
        elseif contains(pokedtext,'left')
            Side = 'Left';
        else
            Side = 'NaN';
        end
        numPresses(t) = numel(pokedIdx);
    else
        Side = 'NaN';
        numPresses(t) = 0;
    end

    
    
    switch TrialType
        case 'FR'
            if isempty(FRside)
                FRside = Side;
            end
        case 'PR'
            if isempty(PRside)
                PRside = Side;
            end
        case 'Miss'
    end
    
    
    
end
TrialList = TrialList';
errorIdx = cellfun(@isempty,TrialList);
TrialList = TrialList(~errorIdx);
numPresses = numPresses';
numPresses = numPresses(~errorIdx);
RewardAmt = str2double(RewardAmt)';
RewardAmt = RewardAmt(~errorIdx);

PRidx = find(get_text(TrialList,'PR'));
FRidx = find(get_text(TrialList,'FR'));
FRrequest = numPresses(FRidx(1));
VolumeRatio = RewardAmt(PRidx(1))/RewardAmt(FRidx(1));

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
    
    PRrequest{y} = PRpresstemp(tempStart:PRlowIdx(tempEnd)-1);
end

% PRpressesunique = unique(sort(PRpresses));
% PRlowval = PRpressesunique
% PRlow = find(PRpresses == 2);
% 
% PRidxx = [0;PRidx];
% prcounter = 0;
% for p = 1:numel(PRlow)
%     PRChk = PRlow(p);
%     for x = PRChk:-1:1
%         prcounter = prcounter+1;
%         if (PRidxx(x) - PRidxx(x-1)) == 1
%         else 
%             break;
%         end
%     end
%     PRtrialsToInceaseByBlock(p) = prcounter;
%     prcounter = 0;
% end
% PRtrialsToIncrease = unique(PRtrialsToInceaseByBlock);
% PRstart = PRlow - PRtrialsToIncrease;
% 
% PRpresstemp = [PRpresses;1];
% PRlower = find(PRpresstemp == 1);
% for pp = 1:numel(PRlow)
%     PRend = find(PRlower > PRlow(pp),1,'first');
%     PRrequest{pp} = PRpresstemp(PRstart(pp):PRlower(PRend)-1);
% end

Info.ID = animalName;
Info.FRSide = FRside;
Info.FRrequest = FRrequest;
Info.PRSide = PRside;
Info.PRtrialsToIncrease = PRTrialsToIncrease;
Info.PRrequest = PRrequest;
Info.VolumeRatio = VolumeRatio;
end