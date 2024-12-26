
function plotPolicyBlock(block)

trials = length(block.HRpress);
Hpress=block.HRpress;
    Hpress(Hpress==0) = NaN;
    Lpress = block.LRpress;
    Lpress(Lpress==0) = NaN;
    Hreward = block.HRreward;
    Hreward(Hreward==0) = NaN;
    Lreward = block.LRreward;
    Lreward(Lreward==0) = NaN;
    block.InitiationTimes=block.InitiationTimes(2:end);
    block.RewardedTimes=block.RewardedTimes(1:end-1);
    ITI=block.InitiationTimes-block.RewardedTimes;
    laser=block.LaserON;
    if nargin<4
        XlimYlim(1,:) = [1 100];
        XlimYlim(2,:) = [0 50];
        XlimYlim(3,:) = [0 3];
        XlimYlim(4,:) = [1 54];
        XlimYlim(5,:)=[1 1000];
    else
    end
    
    figure('Name',[block.File '-']);
    subplot(2,4,[1 2])
    scatter(1:trials,Hpress,150,[0.7 0.5 0.5],'filled','MarkerEdgeColor',[1 1 1],'MarkerFaceAlpha',0.5)
    hold on
    scatter(1:trials,Lpress,150,[0.5 0.5 0.7],'filled','MarkerEdgeColor',[1 1 1],'MarkerFaceAlpha',0.5)
    hold on
    ylim(XlimYlim(2,:))
    xlim(XlimYlim(1,:))
    hold on
    laser= 50*laser;
    plot(laser,'LineWidth',2,'Color',[0.1 0.7 1])
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
    subplot(2,4,[3 4])
    plot(2:trials,ITI,'-','color',[0 0 0],'LineWidth',1.2)

    xlim(XlimYlim(1,:))
    set(gca,'Box','off','FontSize',18,'Box','off','FontName','Dialog','LineWidth',1.2);

    
    subplot(2,4,[5 6])
    scatter(1:trials,Hreward,150,[1 0 0],'filled','MarkerEdgeColor',[1 1 1],'MarkerFaceAlpha',0.5)
    hold on
    scatter(1:trials,Lreward,150,[0 0 1],'filled','MarkerEdgeColor',[1 1 1],'MarkerFaceAlpha',0.5)
    hold on
    Hreward(isnan(Hreward)) = 0;
    Lreward(isnan(Lreward)) = 0;
    HRconfedence = smoothdata(((double(logical(Hreward))-double(logical(Lreward)))+1)/2,'gaussian',6,'omitnan');
    HRconfedence = 2*HRconfedence+1;
    plot(1:trials,HRconfedence,'-','color',[0 0 0],'LineWidth',1.2)
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

    subplot(2,4,[7 8])
    [X_HRrequest,Y_pHR] = estimatePolicy(block.HRrequest,block.HRchoice);
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