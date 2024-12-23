function SEM = estimateSEM(Data)
    % Data  1-D: samples; 2-D: data
    if nargin<2
        BootstrapingSampleNumber = 20;
    else
    end
    if size(Data,1)==1
        Data = Data';
    else
    end
        
    for data = 1:size(Data,2)
        currentData = Data(:,data);
        currentData(isnan(currentData)) = [];
        try
            SEM(1,data) = nanstd(bootstrp(1000,@mean,currentData));
        catch
            SEM(1,data) = 0;
        end
    end
end