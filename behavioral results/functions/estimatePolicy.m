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