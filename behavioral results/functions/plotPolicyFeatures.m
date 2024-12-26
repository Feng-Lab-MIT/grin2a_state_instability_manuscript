function plotPolicyFeatures(PickedBlock,slope1,opto1,slope2,opto2)
    color1 = [0.5294    0.9216    0.8392];
    color2 = [ 0.9216    0.5294    0.8078];
    
    data1 = eval(['PickedBlock.' slope1 '.' opto1 ';']);
    data2 = eval(['PickedBlock.' slope2 '.' opto2 ';']);
    
    X1all = [];
    Y1all = [];
    for block = 1:size(data1,2)
        if contains(slope1,'cont2')
            if max(rem(data1{block}.X_HRrequest(1:end-1),2)==0)==1
                continue
            else
            end
        else
        end
        X1all = [X1all data1{block}.X_HRrequest(1:end-1)]; % (1:end-1) to avoid imprecise end
        Y1all = [Y1all data1{block}.Y_pHR(1:end-1)];
    end
    X1pool = unique(X1all);
    c = 0;
    for request = X1pool
        ypool = Y1all;
        ypool(X1all~=request) = [];
        c = c + 1;
        Y1pool(1,c) = nanmean(smoothdata(ypool,'gaussian',6,'omitnan'));
        E1pool(1,c) = estimateSEM(smoothdata(ypool,'gaussian',6,'omitnan'));
        clear ypool
    end
    
    X2all = [];
    Y2all = [];
    for block = 1:size(data2,2)
        if contains(slope2,'cont2')
            if max(rem(data2{block}.X_HRrequest(1:end-1),2)==0)==1
                continue
            else
            end
        else
        end
        X2all = [X2all data2{block}.X_HRrequest(1:end-1)];
        Y2all = [Y2all data2{block}.Y_pHR(1:end-1)];
    end
    X2pool = unique(X2all);
    c = 0;
    for request = X2pool
        ypool = Y2all;
        ypool(X2all~=request) = [];
        c = c + 1;
        Y2pool(1,c) = nanmean(smoothdata(ypool,'gaussian',6,'omitnan'));
        E2pool(1,c) = estimateSEM(smoothdata(ypool,'gaussian',6,'omitnan'));
        clear ypool
    end
    
    figure('Name',[slope1 ' ' opto1 ' (green) vs ' slope2 ' ' opto2 ' (red)']);
    h = boundedline(X1pool,Y1pool,E1pool,'cmap',color1,'alpha');
    set(h,'LineWidth',1.2)
    hold on
    h = boundedline(X2pool,Y2pool,E2pool,'cmap',color2,'alpha');
    set(h,'LineWidth',1.2)
    hold off
    ylim([0.2 1])
    xlim([1 36])
end