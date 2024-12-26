function Blocks = getBlocks_withITI(Path)
    % E.g. Blocks = arrangeBlocks('C:\Users\Tingting\Google Drive\data_roomF\test');
    addpath(Path)
    Files = dir(fullfile(Path,'*.txt'));
    reverseStr = '';
    for file = 1:length(Files)
        msg = sprintf(['<<< Preprocessing File: ' Files(file).name ' >>>  \n']);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        if contains(Files(file).name,'training') || contains(Files(file).name,'ID') || contains(Files(file).name,'TEST')
            continue
        else
        end
        if exist([Files(file).name(1:end-4) '_info.mat'])
           load([Files(file).name(1:end-4) '_info.mat']);
        else
           Info = [];
        end
        info{file} = getInfoFromTxt([Path '/' Files(file).name],Info);
        if isempty(info{file}.BlockStart) || isempty(info{file}.BlockEnd)
            data{file} = [];
        else
            for block = 1:length(info{file}.BlockEnd)
                startT = info{file}.BlockStart(block);
                endT = info{file}.BlockEnd(block);
                if ~commitPass(info{file})
                    data{file}{block}.File = info{file}.File;
                    data{file}{block}.RewardedTimes = info{file}.RewardedTimes(startT:endT)-nanmin(info{file}.InitiationTimes(startT:endT));
                    data{file}{block}.InitiationTimes = info{file}.InitiationTimes(startT:endT)-nanmin(info{file}.InitiationTimes(startT:endT));
                    data{file}{block}.LRrequest = info{file}.LRrequest;
                    data{file}{block}.HRpattern = info{file}.HRpattern{block};
                    data{file}{block}.VolumeRatio = info{file}.VolumeRatio;
                    data{file}{block}.LRchoice = info{file}.LRchoice(startT:endT);
                    data{file}{block}.HRchoice = info{file}.HRchoice(startT:endT);
                    data{file}{block}.LRpress = info{file}.LRpress(startT:endT);
                    data{file}{block}.HRpress = info{file}.HRpress(startT:endT);
                    data{file}{block}.LRreward = info{file}.LRreward(startT:endT);
                    data{file}{block}.HRreward = info{file}.HRreward(startT:endT);
                    data{file}{block}.LaserON = info{file}.LaserON(startT:endT);
                    data{file}{block}.Manipulation = info{file}.Manipulation;
                    data{file}{block}.SlopeType = info{file}.SlopeType{block};
                    data{file}{block}.HRrequest = info{file}.HRrequest{block};
                    data{file}{block}.PressGoals = info{file}.PressGoals(startT:endT);
                    [data{file}{block}.Policy.X_HRrequest,data{file}{block}.Policy.Y_pHR] = estimatePolicy(data{file}{block}.HRrequest,data{file}{block}.HRchoice);
                    if contains(data{file}{block}.Manipulation,'ds22Q11') || contains(data{file}{block}.Manipulation,'Grin2a') || data{file}{block}.VolumeRatio~=3
                        continue
                    else
                    end
                    ST = data{file}{block}.SlopeType;
                    seperateUnsure = true;
                    if ~seperateUnsure
                        if contains(data{file}{block}.SlopeType,'_unsure')
                            ST(end-6:end) = '';
                        else
                        end
                    else
                    end
                    blockType = ST;
                    if nansum(data{file}{block}.LaserON)~=0
                        laserTrials = data{file}{block}.LaserON;
                        HRchoice = logical(data{file}{block}.HRreward);
                        LRchoice = logical(data{file}{block}.LRreward);
                        rewarded = HRchoice|LRchoice;
                        HRrequest = data{file}{block}.HRrequest;
                        if nansum(laserTrials)==1
                            laserType = [data{file}{block}.Manipulation '_single'];
                        elseif tsnanmode(laserTrials(HRchoice==1))==1 && tsnanmode(laserTrials(LRchoice==1))==0
                            laserType = [data{file}{block}.Manipulation '_HRspesific'];
                        elseif tsnanmode(laserTrials(HRchoice==1))==0 && tsnanmode(laserTrials(LRchoice==1))==1
                            laserType = [data{file}{block}.Manipulation '_LRspesific'];
                        elseif tsnanmode(laserTrials(HRrequest==1&rewarded))==1 && tsnanmode(laserTrials(HRrequest>1&rewarded))==0
                            laserType = [data{file}{block}.Manipulation '_InterBlock'];
                        elseif tsnanmode(laserTrials(HRrequest<19&rewarded))==1 && tsnanmode(laserTrials(HRrequest>19&rewarded))==0
                            laserType = [data{file}{block}.Manipulation '_early'];
                        elseif tsnanmode(laserTrials(HRrequest<19&rewarded))==0 && tsnanmode(laserTrials(HRrequest>19&rewarded))==1
                            laserType = [data{file}{block}.Manipulation '_late'];
                        else
                            laserType = data{file}{block}.Manipulation;
                        end
                    else
                        laserType = 'None';
                    end
                    try
                        eval(['Blocks.' blockType '.' laserType ';'])
                    catch
                        eval(['Blocks.' blockType '.' laserType ' = [];'])
                    end
                    eval(['Blocks.' blockType '.' laserType '{length(Blocks.' blockType '.' laserType ') + 1} = data{file}{block};'])
                else
                end
            end
        end
    end
    function info = getInfoFromTxt(fn,Info)
        behavior = textread(fn,'%s',-1,'delimiter','\t');
        [~,filename] = fileparts(fn);
        filename_parts = split(filename,'_');
        animalNameChk = 2;
        animalID = filename_parts{animalNameChk};
        if ~exist('animalID')
            animalID = ' ';
        end
        get_text = @(x,y) (cellfun(@(c) contains(c,y),x));
        trial_starts = get_text(behavior,'trial available');
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
            if ~isempty(BlockTxt)
                BlockNum(tt) = str2double(BlockTxt{1}(end));
            else
            end
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
                Side = 'Unrecognized';
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
                    Side = 'Unrecognized';
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
            LaserOnChk = get_text(trial_text,'las'); %'laser on'
            if sum(LaserOnChk) == 0
                LaserTrial(tt) = 0;
            elseif sum(LaserOnChk) > 0
                LaserTrial(tt) = 1;
            end
            GoalChk = get_text(trial_text,'goal counter');
            if sum(GoalChk)==0
                PressGoal(tt) = NaN;
            else
                [~,GoalChkIdx] = nanmax(GoalChk);
                if contains(trial_text{GoalChkIdx}(end-1),'=')
                    PressGoal(tt) = str2num(trial_text{GoalChkIdx}(end));
                elseif contains(trial_text{GoalChkIdx}(end-2),'=')
                    PressGoal(tt) = str2num(trial_text{GoalChkIdx}(end-1:end));
                elseif contains(trial_text{GoalChkIdx}(end-3),'=')
                    PressGoal(tt) = str2num(trial_text{GoalChkIdx}(end-2:end));
                else
                    PressGoal(tt) = NaN;
                end
            end
            
            if nanmax(completeChk)==1
                try
                    RewardedTimes(tt) = str2num(trial_text{find(completeChk==1)+1});
                catch
                    RewardedTimes(tt) = NaN;
                end
            else
                RewardedTimes(tt) = NaN;
            end
            
            if nanmax(initChk)==1
                try
                    InitiationTimes(tt) = str2num(trial_text{find(initChk==1)+1});
                catch
                    InitiationTimes(tt) = NaN;
                end
            else
                InitiationTimes(tt) = NaN;
            end
            
            tt = tt+1;
        end
        try
            info.RewardedTimes = (RewardedTimes-nanmin(InitiationTimes))/100;
        catch
            info.RewardedTimes = [];
        end
        try
            info.InitiationTimes = (InitiationTimes-nanmin(InitiationTimes))/100;
        catch
            info.InitiationTimes = [];
        end
        
        TrialList = TrialList';
        if exist('numPresses')
            numPresses = numPresses';
        else
            numPresses = [];
        end
        RewardAmt = str2double(RewardAmt)';
        if exist('BlockNum')
            BlockNum = BlockNum'+1;
            BlockDiff = diff(BlockNum);
            BlockSwitch = find(BlockDiff);
        else
            BlockSwitch = [];
        end
        if exist('rightchoice') && exist('leftchoice')
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
        else
            LRchoice = [];
            HRchoice = [];
            LRpresses = [];
            HRpresses = [];
        end
        LRreward = (strcmp(TrialList,'LR') & ~isnan(RewardAmt));
        HRreward = (strcmp(TrialList,'HR') & ~isnan(RewardAmt));
        PRidx = find(get_text(TrialList,'HR'));
        FRidx = find(get_text(TrialList,'LR'));
        if ~isempty(FRidx) && ~isempty(PRidx)
            FRrequest = numPresses(FRidx(1));
            VolumeRatio = RewardAmt(PRidx(1))/RewardAmt(FRidx(1));
        else
            FRrequest = NaN;
            VolumeRatio = NaN;
        end
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
        if exist('PRTrialsToIncreaseByBlock')
            PRTrialsToIncrease = min(PRTrialsToIncreaseByBlock);
        else
            PRTrialsToIncrease = NaN;
        end
        PRlow = min(PRpresses);
        PRpresstemp = [PRpresses;PRlow];
        PRlowIdx = find(PRlow == PRpresstemp);
        if exist('PRstartIncrease')
            for y = 1:numel(PRstartIncrease)
                tempIncrease = PRstartIncrease(y);
                tempStart = tempIncrease - PRTrialsToIncreaseByBlock(y);
                tempEnd = find(PRlowIdx > tempIncrease,1,'first');
                BlockStartIdx(y) = PRidx(tempStart);
                PRrequest{y} = PRpresstemp(tempStart:PRlowIdx(tempEnd)-1);
            end
        else
            PRrequest = [];
        end
        if exist('BlockStartIdx') && exist('BlockEndIdx')
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
        else
            BlockList = zeros(numel(TrialList),1);
        end
        info.File = filename;
        info.ID = animalID;
        info.LRrequest = FRrequest;
        info.HRSide = PRside;
        info.HRtrialsToIncrease = PRTrialsToIncrease;
        info.HRpattern = PRrequest;
        info.VolumeRatio = VolumeRatio;
        info.LRchoice = LRchoice;
        info.HRchoice = HRchoice;
        info.LRpress = LRpresses;
        info.HRpress = HRpresses;
        info.LRreward = double(LRreward);
        info.HRreward = VolumeRatio.*double(HRreward);
        if exist('PressGoal')
            info.PressGoals = PressGoal';
        else
            info.PressGoals = [];
        end
        if exist('LaserTrial')
            info.LaserON = LaserTrial';
            if nanmax(info.LaserON)==1
                if ~isempty(Info)
                    if isfield(Info,'OptoType') % to add PLtoMD
                        if contains(animalID,'983') || contains(animalID,'984') || contains(animalID,'992')...
                                || contains(animalID,'M0130') || contains(animalID,'M0131')
                            if contains(Info.OptoType,'MD')
                                info.Manipulation = Info.OptoType;
                            else
                                if isfield(Info,'OptoColor') 
                                    if contains(animalID,'M0131')
                                        if str2num(Info.OptoColor(1))==4
                                            info.Manipulation = ['MDto' Info.OptoType];
                                        else
                                            info.Manipulation = Info.OptoType;
                                        end
                                    elseif contains(animalID,'M0130') || contains(animalID,'M0132')
                                        if str2num(Info.OptoColor(1))==4
                                            info.Manipulation = Info.OptoType;
                                        else
                                            info.Manipulation = ['MDto' Info.OptoType];
                                        end
                                    else
                                        info.Manipulation = ['MDto' Info.OptoType];
                                    end
                                else
                                    info.Manipulation = ['MDto' Info.OptoType];
                                end
                            end
                        elseif contains(animalID,'M0122') || contains(animalID,'M0123')
                            if contains(Info.OptoType,'MD')
                                info.Manipulation = [Info.OptoType 'drd2'];
                            else
                                info.Manipulation = ['MDdrd2To' Info.OptoType];
                            end
                        elseif contains(animalID,'M0124') || contains(animalID,'M0125')
                            if contains(Info.OptoType,'MD')
                                info.Manipulation = [Info.OptoType 'grik4'];
                            else
                                info.Manipulation = ['MDgrik4To' Info.OptoType];
                            end
                        else
                            info.Manipulation = Info.OptoType;
                        end
                    else
                        info.Manipulation = 'Unknown';
                    end
                else
                    info.Manipulation = 'Unknown';
                end
            else
                info.Manipulation = 'None';
            end
            if contains(animalID,'M0126') || contains(animalID,'M0127') || contains(animalID,'M0128') || contains(animalID,'M0129')
                info.Manipulation = 'ds22Q11';
            elseif contains(animalID,'TT01') || contains(animalID,'TT08') || contains(animalID,'TT10')
                info.Manipulation = 'Grin2a';
            else
            end
        else
            info.LaserON = [];
            info.Manipulation = 'None';
        end
        info.BlockNumList = BlockList;
        if exist('BlockStart')
            BS = find(BlockStart==1);
        else
            BS = [];
        end
        if exist('BlockStart')
            BE = find(BlockEnd==1);
        else
            BE = [];
        end
        info.BlockStart = [];
        info.BlockEnd = [];
        for test = 1:length(BS)
            if nansum(HRreward(BS(test):BS(test)+PRTrialsToIncrease-1))==PRTrialsToIncrease
                info.BlockStart(test) = BS(test);
                EndTest = BE-BS(test);
                EndTest(EndTest<1) = [];
                if isempty(EndTest)
                    break
                else
                    info.BlockEnd(test) = BS(test)+nanmin(EndTest);
                end
            else
            end
        end

        if length(info.BlockEnd)>0
            slopetypetable{1} = 'cont1'; 
            slopetypetable{2} = 'disc1'; 
            slopetypetable{3} = 'cont2'; 
            slopetypetable{4} = 'disc2'; 
            slopetypetable{5} = 'cont3'; 
            slopetypetable{6} = 'disc3'; 
            slopetypetable{7} = 'x05'; 
            slopetypetable{8} = 'x0'; 
            slopetypetable{9} = 'x4td'; 
            slopetypetable{10} = 'surp_earlyInc_cnt'; 
            slopetypetable{11} = 'surp_earlyDec_cnt'; 
            slopetypetable{12} = 'surp_lateInc_cnt'; 
            slopetypetable{13} = 'surp_lateDec_cnt'; 
            slopetypetable{14} = 'surp_midInc_cnt'; 
            slopetypetable{15} = 'surp_midDec_cnt'; 
            slopetypetable{16} = 'surp_earlyInc_dic'; 
            slopetypetable{17} = 'surp_earlyDec_dic'; 
            slopetypetable{18} = 'surp_lateInc_dic'; 
            slopetypetable{19} = 'surp_lateDec_dic'; 
            slopetypetable{20} = 'surp_midInc_dic'; 
            slopetypetable{21} = 'surp_midDec_dic'; 
            slopetypetable{22} = 'step_MidToHigh'; 
            slopetypetable{23} = 'step_MidToHigh'; 
            slopetypetable{24} = 'surp_midInc_cnt'; 
            slopetypetable{25} = 'step_1to99'; 
            HRsequence{1} = [ones(1,3) 1:1000];
            HRsequence{2} = [ones(1,3) sort([1:2:1000 1:2:1000])];
            HRsequence{3} = [ones(1,3) 1:2:1000];
            HRsequence{4} = [ones(1,3) sort([1:4:1000 1:4:1000])];
            HRsequence{5} = [ones(1,3) 1:4:1000];
            HRsequence{6} = [ones(1,3) sort([1:6:1000 1:6:1000])];
            HRsequence{7} = [ones(1,3) sort([1:1000 1:1000])];
            HRsequence{8} = [ones(1,100)];
            HRsequence{9} = [ones(1,3) sort([1:4:1000 1:4:1000 1:4:1000 1:4:1000])];
            HRsequence{10} = [ones(1,3) 1:1000]; HRsequence{10}(1,15) = 18;
            HRsequence{11} = [ones(1,3) 1:1000]; HRsequence{11}(1,15) = 6;
            HRsequence{12} = [ones(1,3) 1:1000]; HRsequence{12}(1,27) = 30;
            HRsequence{13} = [ones(1,3) 1:1000]; HRsequence{13}(1,27) = 18;
            HRsequence{14} = [ones(1,3) 1:1000]; HRsequence{14}(1,21) = 24;
            HRsequence{15} = [ones(1,3) 1:1000]; HRsequence{15}(1,21) = 12;
            HRsequence{16} = [ones(1,3) sort([1:2:1000 1:2:1000])]; HRsequence{16}(1,15) = 18;
            HRsequence{17} = [ones(1,3) sort([1:2:1000 1:2:1000])]; HRsequence{17}(1,15) = 6;
            HRsequence{18} = [ones(1,3) sort([1:2:1000 1:2:1000])]; HRsequence{18}(1,27) = 30;
            HRsequence{19} = [ones(1,3) sort([1:2:1000 1:2:1000])]; HRsequence{19}(1,27) = 18;
            HRsequence{20} = [ones(1,3) sort([1:2:1000 1:2:1000])]; HRsequence{20}(1,21) = 24;
            HRsequence{21} = [ones(1,3) sort([1:2:1000 1:2:1000])]; HRsequence{21}(1,21) = 12;
            HRsequence{22} = [ones(1,3) 1:17 18.*ones(1,18) 35.*ones(1,100)];
            HRsequence{23} = [ones(1,3) 1:17 18.*ones(1,18) 36.*ones(1,100)];
            HRsequence{24} = [ones(1,3) 1:1000]; HRsequence{24}(1,22) = 36;
            HRsequence{25} = [ones(1,21) 99.*ones(1,100)]; 
            for block = 1:length(info.BlockEnd)
                if nanmax(info.PressGoals(info.BlockStart(block):info.BlockEnd(block)))==99
                    slopetype = 25;
                    info.SlopeType{block} = slopetypetable{slopetype};
                    ForCheck = info.HRpattern{block};
                    shift = 1;
                else
                    CHECK = Inf;
                    maybe = 1;
                    ForCheck = info.HRpattern{block};
                    AlreadyFound = false;
                    for slopetype = 1:size(slopetypetable,2)
                        if ~AlreadyFound
                            shift1 = nansum(double(logical(HRsequence{slopetype}(1:length(ForCheck))-ForCheck')));
                            shift2 = nansum(double(logical(HRsequence{slopetype}(2:length(ForCheck))-ForCheck(1:end-1)')));
                            check = nanmin(shift1,shift2);
                            if check==0
                                info.SlopeType{block} = slopetypetable{slopetype};
                                if shift1<shift2
                                    shift = 1;
                                else
                                    shift = 2;
                                end
                                AlreadyFound = true;
                                continue
                            else
                                if check<CHECK
                                    CHECK = check;
                                    maybe = slopetype;
                                    if shift1<shift2
                                        shift = 1;
                                    else
                                        shift = 2;
                                    end
                                else
                                end
                                if slopetype==size(slopetypetable,1)
                                    info.SlopeType{block} = [slopetypetable{maybe} '_unsure']; % can be within-block
                                else
                                end
                            end
                        else
                        end
                    end
                end
                info.HRrequest{block} = info.HRpress;
                info.HRrequest{block}(info.HRreward==0) = 0;
                info.HRrequest{block} = flip(info.HRrequest{block}(info.BlockStart(block):info.BlockEnd(block)));
                for trial = 1:length(info.HRrequest{block})
                    if trial==1 && info.HRrequest{block}(trial)==0
                        if ~exist('maybe')
                            sequence = HRsequence{slopetype};
                        else
                            sequence = HRsequence{maybe};
                        end
                        if shift==1
                            info.HRrequest{block}(trial) = sequence(length(ForCheck)+1);
                        else
                            info.HRrequest{block}(trial) = sequence(length(ForCheck)+2);
                        end
                    elseif trial>1 && info.HRrequest{block}(trial)==0
                        info.HRrequest{block}(trial) = info.HRrequest{block}(trial-1);
                    else
                    end
                end
                if exist('sequence')
                    info.HRpattern{block} = sequence';
                else
                end
                info.HRrequest{block} = flip(info.HRrequest{block});
                clear sequence
            end
        else
            info.SlopeType = [];
            info.HRrequest = [];
        end
    end
    function pass = commitPass(Info)
        startT = Info.BlockStart(block);
        endT = Info.BlockEnd(block);
        [~,Y_pHR] = estimatePolicy(Info.HRrequest{block},Info.HRchoice(startT:endT));
        if contains(Info.Manipulation,'None') || nansum(info{file}.LaserON(startT:endT))==0
            if nansum(Info.LRreward(endT-5:endT))~=6 || Y_pHR(2)~=1
                pass = true;
                return
            else
                if ~isempty(str2num(Info.SlopeType{block}(end)))
                    if endT-startT+1>100/str2num(Info.SlopeType{block}(end))
                        pass = true;
                        return
                    else
                    end
                else
                end
            end
        else
            if (contains(Info.Manipulation,'MD')&&nansum(Info.LRreward(endT-5:endT))==6)||(contains(Info.Manipulation,'PL')&&sum(Info.LRreward(endT-5:endT))==6)
                pass = true;
                return
            else
            end
        end
        if ~exist('pass')
            pass = false;
        else
        end
    end
    function [X_HRrequest,Y_pHR] = estimatePolicy(HRrequest,HRchoice)
        X_HRrequest = unique(HRrequest);
        c = 0;
        for request = X_HRrequest'
            choice = HRchoice;
            choice(HRrequest~=request) = [];
            c = c + 1;
            Y_pHR(1,c) = nansum(choice)./length(choice);
            if Y_pHR(1,c)==0 && length(choice)<6
                Y_pHR(1,c) = NaN;
            else
            end
            clear choice
        end
        if nanmin(Y_pHR)~=0
            X_HRrequest(end+1) = 150;
            Y_pHR(end+1) = 0;
        else
        end
    end
end