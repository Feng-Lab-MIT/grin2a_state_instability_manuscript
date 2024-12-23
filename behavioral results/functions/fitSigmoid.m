function [k,model,Y_model] = fitSigmoid(Y,X)
    % Input only Y: Y=HRchoice
    % Inputs inlcude Y and X: Y=Y_pHR; X=X_HRrequest
    warning off
    type = 'Normal';
    model = @(k,x)  k(3).*(1-cdf(type,x,k(1),k(2)))+(1-k(3));
    try
        if nargin==1
            if size(Y,2)>size(Y,1)
                Y = Y';
            else
            end
            k = lsqcurvefit(model,[length(Y)/2 3 1],1:length(Y),Y',[6 0.5 0.5],[50 5 1],optimset('Display','off'));
            Y_model = model(k,[1:100]');
        else
            if size(Y,2)>size(Y,1)
                Y = Y';
            else
            end
            if size(X,2)>size(X,1)
                X = X';
            else
            end
            k = lsqcurvefit(model,[18 3 1],X',Y',[-10 0.0001 0.5],[30 5 1],optimset('Display','off'));
            Y_model = model(k,[1:100]');
        end
    catch
        k = nan(1,3);
        Y_model = nan(1,100);
    end
end