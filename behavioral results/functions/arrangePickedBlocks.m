function DATA = arrangePickedBlocks(DATASET)
    % E.g. DATA = arrangePickedBlocks('con1X_02072021')
    load('IDs.mat')
    T = readtable([DATASET '.txt']);
    for a = 1:size(T,1)
        AA = strsplit(T.Var1{a},'-');
        TT(a,1) = str2num(AA{1});
        TT(a,2) = str2num(AA{2});
        TT(a,3) = str2num(AA{3});
    end
    disp(['Arranging ' DATASET '... '])
    reverseStr = '';
    c = 0;
    for bk = 1:size(TT,1)
        msg = sprintf(['  Block %d out of %d  \n'],bk,size(TT,1));
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        animal = TT(bk,1);
        session = TT(bk,2);
        block = TT(bk,3);
        load(['Data_' IDs{animal} '.mat'])
        try
            startT = Data{session}.Trials.BlockStartTrial(block);
            endT = Data{session}.Trials.BlockEndTrial(block);
        catch
            continue
        end
        c = c + 1;
        DATA{c}.HRchoice = Data{session}.Trials.HRchoice(startT:endT);
        DATA{c}.LRchoice = Data{session}.Trials.LRchoice(startT:endT);
        DATA{c}.HRpress = Data{session}.Trials.HRpress(startT:endT);
        DATA{c}.LRpress = Data{session}.Trials.LRpress(startT:endT);
        DATA{c}.HRreward = Data{session}.Trials.HRreward(startT:endT);
        DATA{c}.LRreward = Data{session}.Trials.LRreward(startT:endT);
        [DATA{c}.ExploreOnTrial,DATA{c}.ExploreOffTrial] = getOnOff(DATA{c}.HRreward);
        [DATA{c}.X_HRrequest,DATA{c}.Y_pHR] = estimatePolicyV4(animal,session,block);
    end
end