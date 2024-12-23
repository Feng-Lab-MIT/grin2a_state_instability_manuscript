function plotSession(Session)


trials=length(Session.HRchoice);
Hpress=Session.HRpress;
Hpress(Hpress==0) = NaN;
Lpress = Session.LRpress;
LGoals =Lpress;
LGoals(LGoals<=6)=6;
Lpress(Lpress==0) = NaN;
Hreward = Session.HRreward;
Hreward(Hreward==0) = NaN;
Lreward = Session.LRreward;
Lreward(Lreward==0) = NaN;
Lchoice=Session.LRchoice;
Hchoice=Session.HRchoice;
PressGoals=Session.PressGoals;
PressGoals(Hchoice==0)=0;

for j=1:length(PressGoals)
    PressGoals(1)=1;
    if PressGoals(j)==0
        PressGoals(j)=PressGoals(j-1);
    end
end
PressGoals(isnan(PressGoals))=0;
Hvalue=3./PressGoals;
Lvalue=1./LGoals;
Session.InitiationTimes=Session.InitiationTimes(2:end);
Session.RewardedTimes=Session.RewardedTimes(1:end-1);
ITI=Session.InitiationTimes-Session.RewardedTimes;
laser=Session.LaserON;

figure('Name',[Session.File '-']);
subplot(4,1,1)
plot(1:trials,log(Hvalue),'LineWidth',2,'Color',[0.8 0.5 0.5])
hold on;
plot(1:trials,log(Lvalue),'LineWidth',2,'Color',[0.5,0.5,0.8])
hold off;
set(gca,'Box','off','FontSize',10,'Box','off','FontName','Dialog','LineWidth',1.0);
xlim([1,trials])
ylabel('log(Value)')

subplot(4,1,2)
scatter(1:trials,Hreward,20,[0.8 0.5 0.5],'filled')
hold on
scatter(1:trials,Lreward,20,[0.5 0.5 0.8],'filled')
hold on
    Hreward(isnan(Hreward)) = 0;
    Lreward(isnan(Lreward)) = 0;
    HRconfedence = smoothdata(((double(logical(Hreward))-double(logical(Lreward)))+1)/2,'gaussian',6,'omitnan');
    HRconfedence = 2*HRconfedence+1;
    plot(1:trials,HRconfedence,'-','color',[0 0 0],'LineWidth',1.2)
    
    hold on
    plot(1:trials,laser+1,'LineWidth',2,'Color',[0,0,1]);
    xlim([1,trials])
    
 set(gca,'Box','off','FontSize',10,'Box','off','FontName','Dialog','LineWidth',1.0);   
    
   
    subplot(4,1,3)
plot(1:trials,PressGoals./3,'LineWidth',2,'Color',[0.8 0.5 0.5])
hold on;
plot(1:trials,LGoals,'LineWidth',2,'Color',[0.5,0.5,0.8])
hold off;
set(gca,'Box','off','FontSize',10,'Box','off','FontName','Dialog','LineWidth',1.0);
xlim([1,trials])
ylabel('cost/benefit')

subplot(4,1,4)
scatter(1:trials,Hpress,20,[0.8 0.5 0.5],'filled')
    hold on
scatter(1:trials,Lpress,20,[0.5 0.5 0.8],'filled')
xlim([1,trials])
set(gca,'Box','off','FontSize',10,'Box','off','FontName','Dialog','LineWidth',1.0);
hold off

xlabel('trials')
end
