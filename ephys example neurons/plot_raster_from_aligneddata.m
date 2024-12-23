function [counts,edges] = plot_raster_from_aligneddata(aligneddata,xleft,xright,sortedIdx,Binning,offset)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%20220507 add linewidth to 1.1
TrialRange=[xleft xright];
edges=TrialRange(1):Binning:TrialRange(2);
counts=zeros(length(sortedIdx),length(edges)-1);


leftend=xleft;
rightend=xright;


for row=1:length(sortedIdx)

        aligndata=aligneddata{sortedIdx(row)};

        if isempty(aligndata)==0
            if length(aligndata)~=2
                plot([aligndata aligndata],[row-1+offset row+offset],'k','LineWidth',1.2);
                hold on;
            else
                plot([aligndata(1) aligndata(1)],[row-1+offset row+offset],'k','LineWidth',1.2);
                plot([aligndata(2) aligndata(2)],[row-1+offset row+offset],'k','LineWidth',1.2);
                hold on;
            end
        end
        counts(row,:)=histcounts(aligndata(aligndata<TrialRange(2) & aligndata>TrialRange(1)),edges)./Binning;    
    
end
%flip y axis
set(gca, 'YDir','reverse')
xlim(TrialRange);
%countssmooth=smoothdata(counts,'gaussian',Smoothing);
end

