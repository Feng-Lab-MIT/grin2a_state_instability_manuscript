function [X_HRrequest,Y_pHR] = estimatePolicyV4(animal,session,block)
   try
        load('IDs.mat')
        load(['Data_' IDs{animal} '.mat'])
        startT = Data{session}.Trials.BlockStartTrial(block);
        endT = Data{session}.Trials.BlockEndTrial(block);
        HRrequest = Data{session}.Trials.HRrequest{block}';
        HRchoice = Data{session}.Trials.HRchoice(startT:endT);
        
        % block-spesific debug
        if animal==3 && session==141 && block==3
            HRrequest(HRrequest==74) = 37;
        else
        end
        % block-spesific debug
        
        while HRrequest(end)==1
            HRrequest = HRrequest(1:end-1);
            if length(unique(HRrequest))==1
                X_HRrequest = NaN;
                Y_pHR = NaN;
                return
            else
            end
            endT = endT - 1;
            HRchoice = Data{session}.Trials.HRchoice(startT:endT);
        end
        X_HRrequest = unique(HRrequest);
        c = 0;
        for request = X_HRrequest
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
        if mode(X_HRrequest(2:end)-X_HRrequest(1:end-1))==1
            X_HRrequest(end) = X_HRrequest(end-1)+1;
        elseif mode(X_HRrequest(2:end)-X_HRrequest(1:end-1))==2
            X_HRrequest(end) = X_HRrequest(end-1)+2;
        else
            X_HRrequest(end) = X_HRrequest(end-1)+1;
        end
        
        
        
        
        
        if nanmin(Y_pHR)~=0
            X_HRrequest(end+1) = 150;
            Y_pHR(end+1) = 0;
        else
            if mode(X_HRrequest(2:end)-X_HRrequest(1:end-1))==1
                addHRrequest = max(X_HRrequest):150;
                addpHR = zeros(1,length(addHRrequest));
                X_HRrequest = [X_HRrequest addHRrequest];
                Y_pHR = [Y_pHR addpHR];
            elseif mode(X_HRrequest(2:end)-X_HRrequest(1:end-1))==2
                addHRrequest = max(X_HRrequest):2:150;
                addpHR = zeros(1,length(addHRrequest));
                X_HRrequest = [X_HRrequest addHRrequest];
                Y_pHR = [Y_pHR addpHR];
            else
                addHRrequest = max(X_HRrequest):150;
                addpHR = zeros(1,length(addHRrequest));
                X_HRrequest = [X_HRrequest addHRrequest];
                Y_pHR = [Y_pHR addpHR];
            end
        end
    catch
        X_HRrequest = NaN;
        Y_pHR = NaN;
    end
end