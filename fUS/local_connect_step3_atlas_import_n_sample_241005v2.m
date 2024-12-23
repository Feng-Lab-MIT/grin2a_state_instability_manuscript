%for resample



VMB=niftiread('\\Fenglab03\Yiyun\20230410_fUS_seedmap\regions\313_MB.nii\313_MB.nii');




%%
vqMB=Resample(VMB);



%%
save('\\Fenglab03\Yiyun\20230410_fUS_seedmap\region_atlas_newly_sample_v5_241005MBresam_newreg.mat','vqMB');



%%

function vqACC = Resample(Vx)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
X=[1.35:-0.0025:-0.0025*1319]; %posterior(more) anterior 1.4 is the X(0) v4:1.2 v5/v6: 1.35  v7:1.4   v8:1.45  v9:1.55
Y=[1.3:-0.0025:-0.0025*800]; %up(more) down 1.4 is the Y(0) %v2:1.2, v3:1.1 v4:1.2 v5/v6:1.3  v7:1.35 v8: 1.4   v9:1.5
Z=[-569*0.0025:0.0025:570*0.0025];


%[xq,yq,zq] = meshgrid(-0.6:0.02:0.6,-1.4:0.05:1.6,-0.8:0.02:1.4); % you may need to adjust here
% X=[1.4:-0.0025:-0.0025*1319];
% Y=[1.4:-0.0025:-0.0025*800];
% Z=[-569*0.0025:0.0025:570*0.0025];
% F = scatteredInterpolant(transformeddata(:,1),transformeddata(:,2),transformeddata(:,3),transformeddata(:,4));
% [xq,yq,zq] = meshgrid(-0.6:0.02:0.6,-1.4:0.05:1.6,-0.8:0.02:1.4); % -0.6
% to 0.6=> lateral, -1.4 to 1.4 anterior posterior
% vq = F(xq,yq,zq);

[~,zi]=min(abs(Z-(-0.6)));
zintval=0.03/0.0025;
[~,zf]=min(abs(Z-(0.6)));

[~,xi]=min(abs(X-(-1.2)));
xintval=0.03/0.0025;
[~,xf]=min(abs(X-(1.2)));

[~,yi]=min(abs(Y-(-0.8)));
yintval=0.03/0.0025;
[~,yf]=min(abs(Y-(1.4)));

VACCgrid=zeros((length(xi:xintval:xf))*(length(yi:yintval:yf))*(length(zi:zintval:zf)),4);


xgrid=xf:xintval:xi;
ygrid=yf:yintval:size(Vx,2);
zgrid=zi:zintval:zf;

m=1;

for i=1:length(xgrid)
    for j=1:length(ygrid)
        for k=1:length(zgrid)
            VACCgrid(m,:)=[Z(zgrid(k)),X(xgrid(i)),Y(ygrid(j)),double(Vx(xgrid(i),ygrid(j),zgrid(k)))];
            m=m+1;
        end
    end
end

[xq,yq,zq] = meshgrid(-0.6:0.03:0.6,-1.2:0.03:1.2,-0.8:0.03:1.4); 

Facc = scatteredInterpolant(VACCgrid(:,1),VACCgrid(:,2),VACCgrid(:,3),VACCgrid(:,4));
vqACC = Facc(xq,yq,zq);


end



