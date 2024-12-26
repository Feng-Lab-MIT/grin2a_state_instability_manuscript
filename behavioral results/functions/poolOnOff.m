function [Onset,Offset]=poolOnOff(Blocks,slope,opto)

data=eval(['Blocks.' slope '.' opto]);
Onset=[];
Offset=[];
for i=1:length(data)
[Onset(i),Offset(i)]=getOnset(data{i}.LRpress);

end

function [Onset, Offset] = getOnset(LRpress)
     onsetseq=find(LRpress==6);
    if length(onsetseq)>=2
    Onset=onsetseq(1);
    else
        Onset=0;
    end
    offsetseq=findConSame(LRpress,6,6);
    Offset=offsetseq(1);
end

function t=findConSame(array,val,n)   %find n consecutive same values in array
%for example, find trial number that animal first got rewarded 6 times consecutively on FR side
% array=[1 0 1 0 0 0 1 1 1 0 0 0 0 0 1 1 1 1 1 0 1 0 1 1 1 1 1 1 1 1 0 0 1 0], val=1, n=6

sub_array={};
T=[];

if length(array)>n
   
for j=1:length(array)-n+1
    sub_array{j}=array(j:j+n-1);
    if unique(sub_array{j})==val%|unique(sub_array{j})==[val,val+2]%&unique(sub_array{j+1})==val+2;
        T(end+1)=j;
    end
end
%for j=1:length(array)-n+1
%    sub_array{j}=array(j:j+n-1);
 %   if unique(sub_array{j})==val
  %      T(end+1)=j;
   % end
%end
if ~isempty(T)
     %t=min(T);
     t=T;
else t=0;
end
else t=0;
end
end
end


    
