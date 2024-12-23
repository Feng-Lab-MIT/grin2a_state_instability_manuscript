function s=ste(data)
%%%%%%%for computing stardard error
[m,n]=size(data);
s=std(data);
s=s/m^0.5;