function Optimality = estimateOptimality(Data,session,block)
    Hpress = Data{session}.HRpress(Data{session}.BlockStart(block):Data{session}.BlockEnd(block));
    Hpress(Hpress==0) = NaN;
    Lpress = Data{session}.LRpress(Data{session}.BlockStart(block):Data{session}.BlockEnd(block));
    Lpress(Lpress==0) = NaN;
    Hreward = Data{session}.HRreward(Data{session}.BlockStart(block):Data{session}.BlockEnd(block));
    Hreward(Hreward==0) = NaN;
    Lreward = Data{session}.LRreward(Data{session}.BlockStart(block):Data{session}.BlockEnd(block));
    Lreward(Lreward==0) = NaN;
    Hpress(isnan(Hpress)) = 0;
    Lpress(isnan(Lpress)) = 0;
    Hreward(isnan(Hreward)) = 0;
    Lreward(isnan(Lreward)) = 0;
    Hpattern = Data{session}.HRpattern{block};
    AllPress_best = Hpattern((nanmax(Hreward)./Hpattern)>(nanmax(Lreward)./Data{session}.LRrequest));
    if nanmax(Data{session}.HRrequest{block})<nanmax(AllPress_best)
        [~,AtItsMax] = nanmax(AllPress_best==nanmax(Data{session}.HRrequest{block}));
        AllPress_best = AllPress_best(1:AtItsMax);
        NeedReduce = 1;
    else
        NeedReduce = 0;
    end
    AllReward_best = nanmax(Hreward).*ones(length(Hpattern((nanmax(Hreward)./Hpattern)>(nanmax(Lreward)./Data{session}.LRrequest))),1); 
    if NeedReduce==1
        AllReward_best = AllReward_best(1:AtItsMax);
    else
    end
    AllPress = Hpress+Lpress;
    AllPress(isnan(AllPress)) = 0;
    AllReward = Hreward+Lreward;
    AllReward(isnan(AllReward)) = 0;
    MeanValue = AllReward./AllPress;
    MeanValue(isnan(MeanValue)) = 0;
    Optimality = nanmean(MeanValue)./nanmean(AllReward_best./AllPress_best);
    Optimality(isinf(Optimality)|Optimality<0|Optimality>1) = NaN;
end