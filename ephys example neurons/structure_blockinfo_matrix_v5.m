function [Blockinfo] = structure_blockinfo_matrix_v5(rasterfilepath,BlockN)
%UNTITLED3 Summary of this function goes here
%   updated 20211007 to modify the order of the block <= better version,
%   fixed the errors in v3

g=load(rasterfilepath);

HRLRchoice=(g.RasterData.SessionInfo.Trials.HRchoice(g.RasterData.SessionInfo.Trials.BlockNumList==BlockN)-g.RasterData.SessionInfo.Trials.LRchoice(g.RasterData.SessionInfo.Trials.BlockNumList==BlockN));
HRLRPressN=(g.RasterData.SessionInfo.Trials.HRpresses(g.RasterData.SessionInfo.Trials.BlockNumList==BlockN)+g.RasterData.SessionInfo.Trials.LRpresses(g.RasterData.SessionInfo.Trials.BlockNumList==BlockN));
HRLRrewardif=(g.RasterData.SessionInfo.Trials.HRreward(g.RasterData.SessionInfo.Trials.BlockNumList==BlockN)-g.RasterData.SessionInfo.Trials.LRreward(g.RasterData.SessionInfo.Trials.BlockNumList==BlockN));


choicetrl=((g.RasterData.SessionInfo.Trials.HRchoice(g.RasterData.SessionInfo.Trials.BlockNumList==BlockN)-g.RasterData.SessionInfo.Trials.LRchoice(g.RasterData.SessionInfo.Trials.BlockNumList==BlockN))~=0);
presstrl=((g.RasterData.SessionInfo.Trials.HRpresses(g.RasterData.SessionInfo.Trials.BlockNumList==BlockN)+g.RasterData.SessionInfo.Trials.LRpresses(g.RasterData.SessionInfo.Trials.BlockNumList==BlockN))~=0);
rewardtrl=((g.RasterData.SessionInfo.Trials.HRreward(g.RasterData.SessionInfo.Trials.BlockNumList==BlockN)-g.RasterData.SessionInfo.Trials.LRreward(g.RasterData.SessionInfo.Trials.BlockNumList==BlockN))~=0);
HRLRchoice=HRLRchoice(rewardtrl & (choicetrl & presstrl));
HRLRPressN=HRLRPressN(rewardtrl & (choicetrl & presstrl));
HRLRrewardif=HRLRrewardif(rewardtrl & (choicetrl & presstrl));
    

    
   
%assert(length(HRLRchoice)==length(HRLRPressN),'length of HRLRchoice is not the same as HRLRpressN');
%assert(length(HRLRchoice)==length(HRLRrewardif),'length of HRLRchoice is not the same as HRLRreward');

Blockinfo=[HRLRchoice,HRLRPressN,HRLRrewardif];
%%

end

