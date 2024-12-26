function PlotBarwithError(x,data1,data2)

if x==2
data1_mean=mean(data1);
data1_sem=estimateSEM(data1);
data2_mean=mean(data2);
data2_sem=estimateSEM(data2);

bar(1:2, [data1_mean, data2_mean]);
hold on
er = errorbar(1:2,[data1_mean, data2_mean],[data1_sem, data2_sem]);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
end


% plot individual data in the plot
hold on;
temp1=ones(size(data1,1),1);
scatter(temp1,data1);

hold on;
temp2=2*ones(size(data2,1),1);
scatter(temp2,data2);

end

