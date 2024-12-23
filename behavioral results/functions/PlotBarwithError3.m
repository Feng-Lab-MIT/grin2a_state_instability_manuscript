function PlotBarwithError3(x,data1,data2,data3)
if x==3
data1_mean=mean(data1);
data1_sem=estimateSEM(data1);
data2_mean=mean(data2);
data2_sem=estimateSEM(data2);

data3_mean=mean(data3);
data3_sem=estimateSEM(data3);

bar(1:3, [data1_mean, data2_mean,data3_mean]);
hold on
er = errorbar(1:3,[data1_mean, data2_mean,data3_mean],[data1_sem,data2_sem,data3_sem]);
hold off
er.Color = [0 0 0];                            
er.LineStyle = 'none';     
end
    
end