function plotBlockPolicy(animal,session,block,XlimYlim)
    load('IDs.mat')
    load(['Data_' IDs{animal} '.mat'])
    startT = Data{session}.Trials.BlockStartTrial(block);
    endT = Data{session}.Trials.BlockEndTrial(block);
            
    trials = startT:endT;
    Hpress = Data{session}.Trials.HRpress(startT:endT);
    Hpress(Hpress==0) = NaN;
    Lpress = Data{session}.Trials.LRpress(startT:endT);
    Lpress(Lpress==0) = NaN;
    Hreward = Data{session}.Trials.HRreward(startT:endT);
    Hreward(Hreward==0) = NaN;
    Lreward = Data{session}.Trials.LRreward(startT:endT);
    Lreward(Lreward==0) = NaN;
    if nargin<4
        XlimYlim(1,:) = [1 100];
        XlimYlim(2,:) = [0 50];
        XlimYlim(3,:) = [0 3];
        XlimYlim(4,:) = [1 54];
    else
    end
    
    figure('Name',[num2str(animal) '-' num2str(session) '-' num2str(block)]);
    subplot(2,4,[1 2])
    scatter(1:length(trials),Hpress,150,[0.7 0.5 0.5],'filled','MarkerEdgeColor',[1 1 1],'MarkerFaceAlpha',0.5)
    hold on
    scatter(1:length(trials),Lpress,150,[0.5 0.5 0.7],'filled','MarkerEdgeColor',[1 1 1],'MarkerFaceAlpha',0.5)
    hold on
    ylim(XlimYlim(2,:))
    xlim(XlimYlim(1,:))
    hold off
    set(gca,'Box','off','FontSize',18,'Box','off','FontName','Dialog','LineWidth',1.2);
    if XlimYlim(1,1)==1 && XlimYlim(1,2)==100
        set(gca,'XTick',[1 10:10:100],'XTickLabel',[1 10:10:100]);
    else
    end
    if XlimYlim(2,1)==0 && XlimYlim(2,2)==50
        set(gca,'YTick',0:10:50,'YTickLabel',0:10:50);
    else
    end

    subplot(2,4,[5 6])
    scatter(1:length(trials),Hreward,150,[1 0 0],'filled','MarkerEdgeColor',[1 1 1],'MarkerFaceAlpha',0.5)
    hold on
    scatter(1:length(trials),Lreward,150,[0 0 1],'filled','MarkerEdgeColor',[1 1 1],'MarkerFaceAlpha',0.5)
    hold on
    Hreward(isnan(Hreward)) = 0;
    Lreward(isnan(Lreward)) = 0;
    HRconfedence = smoothdata(((double(logical(Hreward))-double(logical(Lreward)))+1)/2,'gaussian',6,'omitnan');
    HRconfedence = 2*HRconfedence+1;
    plot(1:length(trials),HRconfedence,'-','color',[0 0 0],'LineWidth',1.2)
    hold off
    ylim(XlimYlim(3,:))
    xlim(XlimYlim(1,:))
    set(gca,'Box','off','FontSize',18,'Box','off','FontName','Dialog','LineWidth',1.2);
    if XlimYlim(1,1)==1 && XlimYlim(1,2)==100
        set(gca,'XTick',[1 10:10:100],'XTickLabel',[1 10:10:100]);
    else
    end
    if XlimYlim(3,1)==0 && XlimYlim(3,2)==3
        set(gca,'YTick',0:3,'YTickLabel',0:3);
    else
    end

    subplot(2,4,[3 4 7 8])
    [X_HRrequest,Y_pHR] = estimatePolicyV4(animal,session,block);
    for request = 1:length(X_HRrequest)
        scatter(X_HRrequest(request),Y_pHR(request),150,[Y_pHR(request) 0 (1-Y_pHR(request))],'filled','MarkerEdgeColor',[1 1 1],'MarkerFaceAlpha',0.5)
        hold on
    end
    ylim([0 1])
    xlim(XlimYlim(4,:))
    set(gca,'Box','off','FontSize',18,'Box','off','FontName','Dialog','LineWidth',1.2);
    if XlimYlim(4,1)==1 && XlimYlim(4,2)==54
        set(gca,'XTick',[0 1 6:6:54],'XTickLabel',[0 1 6:6:54]);
    else
    end

    set(gcf,'color','w','Position',[230 513 1158 420]); 
end