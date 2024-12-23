function [Onset,Offset] = getOnset(Data,session,block)


    LRchoice = Data{session}.LRchoice(Data{session}.BlockStart(block):Data{session}.BlockEnd(block));
%     HRchoice = Data{session}.HRchoice(Data{session}.BlockStart(block):Data{session}.BlockEnd(block));
    
    onsetseq=find(LRchoice==1);
    if length(onsetseq)>=2
    Onset=onsetseq(1);
    else
        Onset=0;
    end
    offsetseq=findConSame(LRchoice,1,6);
    Offset=offsetseq(1);
end

