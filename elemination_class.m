clc
clear

s2Wav=[443 490 560 665 705 740 783 842 865 945 1610 2190];
data1=readcell("G:\Trainee\pragati\Copy of modified data_nal.xlsx");
ref1=cell2mat(data1(2:end,5:16));
class1=data1(2:end,25);

ref00=ref1./10000;
Wclass00=class1;
uc=unique(Wclass00);
nonuse=strcmp(Wclass00,uc(5));
ref00(nonuse,:)=[];
Wclass00(nonuse,:)=[];


data1=readcell("G:\Trainee\pragati\Nal_100_DataPoints.xlsx");
ref1=cell2mat(data1(2:end,5:16));
class1=data1(2:end,25);
ref01=ref1./10000;
Wclass01=class1;

ref=[ref00;ref01];
Wclass=[Wclass00;Wclass01];



uc=unique(Wclass);
figure;axes1 = axes(gcf);set(axes1,'FontSize',18);hold on;box on
hold on
for i=1:length(uc)
    plot(nanmean(ref(strcmp(Wclass,uc(i)),:)),LineWidth=2)
end
legend(uc,Location="best")
xlabel('Sentinel-2 Band #');xlim([1 12])
ylabel('Reflectance')


colorpallet=[118 171 41;216 82 24;236 176 31;125 46 141;0 113 188]/255;
figure;axes1 = axes(gcf);set(axes1,'FontSize',18);hold on;box on
hold on
for i=1:length(uc)
    errorbar(nanmean(ref(strcmp(Wclass,uc(i)),:)),nanstd(ref(strcmp(Wclass,uc(i)),:)),LineWidth=2,Color=colorpallet(i,:))
end
legend(uc,Location="best")
xlabel('Sentinel-2 Band #');xlim([1 12])
ylabel('Reflectance')

figure;axes1 = axes(gcf);set(axes1,'FontSize',18);hold on;box on
hold on
for i=1:length(uc)
    plot((ref(strcmp(Wclass,uc(i)),:)'),'-.',LineWidth=0.25,Color=colorpallet(i,:))
    plot(nanmean(ref(strcmp(Wclass,uc(i)),:)),LineWidth=2,Color=colorpallet(i,:))
end
legend(uc,Location="best")
xlabel('Sentinel-2 Band #');xlim([1 12])
ylabel('Reflectance')




maj5=sum([strcmp(Wclass,'Algae') strcmp(Wclass,'Emergent') strcmp(Wclass,'Submerged') strcmp(Wclass,'Water') strcmp(Wclass,'Land')],2);
maj5in=find(maj5==1);





NDVI=(ref(:,8)-ref(:,4))./(ref(:,8)+ref(:,4)); % Normalized Difference Vegetation Index
MNDWI1=(ref(:,3)-ref(:,11))./(ref(:,3)+ref(:,11)); % Modified Normalized Difference Water Index using SWIR band 1
NDWI=(ref(:,3)-ref(:,8))./(ref(:,3)+ref(:,8)); %  Normalized Difference Water Index
NDTI=(ref(:,4)-ref(:,3))./(ref(:,3)+ref(:,4)); %  Normalized Difference Turbidity Index
NDMI1=(ref(:,8)-ref(:,11))./(ref(:,8)+ref(:,11)); %  Normalized Difference Moisture Index using SWIR band 1
EVI = 2.5.*(ref(:,8)-ref(:,4))./(ref(:,8)+6.*ref(:,4)-7.5.*ref(:,2)); % Enhenced Vegetation Index
FAI =  ref(:,8) - (ref(:,4) + (ref(:,11)-ref(:,4))*0.1873); % Floating Algae Index
AFAI = ref(:,8) - ref(:,4) + (ref(:,11)-ref(:,4))*0.5; % Alternative Floating Algae Index
ABDI = (ref(:,6)-ref(:,4))-(ref(:,8)-ref(:,4))*0.375 - (ref(:,4)-0.5*ref(:,3)); % Algal Bloom Detection Index
TMI = (50*(ref(:,3)-ref(:,4)).*ref(:,8))./ref(:,4); % Three band Macro-Algae detection Index
SABI = (ref(:,8)-ref(:,4))./(ref(:,2)+ref(:,3)); % Surface Algal Bloom Index
VBFAH = ((ref(:,8)-ref(:,3))+(ref(:,3)-ref(:,4)).*(ref(:,8)-ref(:,3)))./(2*ref(:,8)-ref(:,4)-ref(:,3)); % vertual baseline floating macro algae height
SA = (ref(:,4)-ref(:,8))/.177; % Slope Algorithm
MAI = ref(:,3) + ref(:,4) - (ref(:,11)+(ref(:,2)-ref(:,11))*0.9375); % Macro-Algae Index
SUI = (ref(:,3)-(ref(:,4)+(ref(:,2)-ref(:,4))*1.5))./ref(:,3); % Sargassum Ulva prolifera Index
NDCI = (ref(:,5)-ref(:,4))./(ref(:,5)+ref(:,4)); % Normalized Difference Chlorophyll Index
MOSES3B = (1./ref(:,4)+1./ref(:,5)).*ref(:,6); % Three band Chlorophyll Index using RE1
S23BDA = (1./ref(:,4)+1./ref(:,5)).*ref(:,8); % Three band Chlorophyll Index using NIR
B3B2 = (ref(:,3)-ref(:,2))./(ref(:,3)+ref(:,2)); % Normalized difference Blue green Index
GNDVI = (ref(:,8)-ref(:,3))./(ref(:,8)+ref(:,3)); % Green Normalized Difference Vegetation Index
EVSI = (ref(:,4)-ref(:,11))./(ref(:,4)+ref(:,11)); % Enhanced Vegetation Index
CMI =  ref(:,3) - (ref(:,2) - (ref(:,11)-ref(:,2))*0.0625); % Cyanobacteria and Macrophytes Index
NDAVI=(ref(:,8)-ref(:,2))./(ref(:,8)+ref(:,2)); % Normalized Difference Aquatic Vegetation Index
WAVI = (ref(:,8)-ref(:,2))./(ref(:,8)+ref(:,2)+0.5); % Water Adjusted Vegetation Index



Windices=[NDVI,MNDWI1,NDWI,NDTI,NDMI1,...
    EVI,FAI,AFAI,ABDI,TMI,...
    SABI,VBFAH,SA,MAI,SUI,...
    NDCI,MOSES3B,S23BDA,B3B2,GNDVI,...
    EVSI,CMI,NDAVI,WAVI];

cl25=Wclass(maj5in);
in25= Windices(maj5in,:);

cl25(6)=[];
in25(6,:)=[];

cl25(204)=[];
in25(204,:)=[];

clearvars -except cl25 in25

[a,~,c]=unique(cl25);
figure;
x=-5:0.001:5;
for i=1:24 
    subplot(6,4,i);
    box on
    hold on
    for j=1:5
    mld=fitdist(in25(c==j,i),"Kernel");
    y=mld.pdf(x);
    p=area(x,y);
%     p(1).EdgeColor=p(1).FaceColor;
    p(1).FaceAlpha=0.5;
    end
    xmin=max(-5,quantile(in25(:,i),0.01));
    xmax=min(5,quantile(in25(:,i),0.99));
    xlim([xmin xmax])
    xlabel('Index Value');
    ylabel('Density')
    hold off
end

%%
clearvars -except cl25 in25
cl5=cl25(:,1);
in5=in25(:,:);
ca=corrcoef(in5);

for i = 1 : 5
    [a,~,c]=unique(cl5);
    c(c==i)=10;
    c(c~=10)=11;
    for j=1:size(in25,2)
        mld = fitctree(in5(:,j),c,'MaxNumSplits',1);
        prd = mld.predict(in5(:,j));
        OA(i,j) = 100*sum(c==prd)/length(c);
    end
end
[a1,b1]=find(OA==max(OA(:)));
a1=a1(1);b1=b1(1);
disp(max(OA(:)))
disp(a(a1))
disp(b1)
clearvars OA
% gmm=fitgmdist(in5(:,b1),2);
% P = posterior(gmm,(-1:.001:1)');
% [~,w]=(min(abs(P(:,1)-0.5)));
% ThresholdG0=-1+w*0.001
[a,~,c]=unique(cl5);
c(c==a1)=10;
c(c~=10)=11;
mld = fitctree(in5(:,b1),c,'MaxNumSplits',1);
view(mld,'Mode','graph')
prd = mld.predict(in5(:,b1));
OA= 100*sum(c==prd)/length(c);


in5(strcmp(cl5,a(a1)),:)=[];
cl5(strcmp(cl5,a(a1))==1)=[];
for i = 1 : 4
    [a,~,c]=unique(cl5);
    c(c==i)=10;
    c(c~=10)=11;
    for j=1:size(in25,2)
        mld = fitctree(in5(:,j),c,'MaxNumSplits',1);
        prd = mld.predict(in5(:,j));
        OA(i,j) = 100*sum(c==prd)/length(c);
    end
end
[a1,b1]=find(OA==max(OA(:)));
a1=a1(1);b1=b1(1);
disp(max(OA(:)))
disp(a(a1))
disp(b1)
clearvars OA
% gmm=fitgmdist(in5(:,b1),2);
% P = posterior(gmm,(-1:.01:1)');
% [~,w]=(min(abs(P(:,1)-0.5)));
% ThresholdG0=-1+w*0.01
[a,~,c]=unique(cl5);
c(c==a1)=10;
c(c~=10)=11;
mld = fitctree(in5(:,b1),c,'MaxNumSplits',1);
view(mld,'Mode','graph')
prd = mld.predict(in5(:,b1));
OA= 100*sum(c==prd)/length(c);


in5(strcmp(cl5,a(a1)),:)=[];
cl5(strcmp(cl5,a(a1))==1)=[];
for i = 1 : 3
    [a,~,c]=unique(cl5);
    c(c==i)=10;
    c(c~=10)=11;
    for j=1:size(in25,2)
        mld = fitctree(in5(:,j),c,'MaxNumSplits',1);
        prd = mld.predict(in5(:,j));
        OA(i,j) = 100*sum(c==prd)/length(c);
    end
end
[a1,b1]=find(OA==max(OA(:)));
a1=a1(1);b1=b1(1);
disp(max(OA(:)))
disp(a(a1))
disp(b1)
clearvars OA
gmm=fitgmdist(in5(:,b1),2);
P = posterior(gmm,(-1:.01:1)');
[~,w]=(min(abs(P(:,1)-0.5)));
ThresholdG0=-1+w*0.01;
[a,~,c]=unique(cl5);
c(c==a1)=10;
c(c~=10)=11;
mld = fitctree(in5(:,b1),c,'MaxNumSplits',1);
view(mld,'Mode','graph')
prd = mld.predict(in5(:,b1));
OA= 100*sum(c==prd)/length(c);


in5(strcmp(cl5,a(a1(1))),:)=[];
cl5(strcmp(cl5,a(a1(1)))==1)=[];
[a,~,c]=unique(cl5);
c(c==1)=10;
c(c~=10)=11;
for j=1:size(in25,2)
    mld = fitctree(in5(:,j),c,'MaxNumSplits',1);
    prd = mld.predict(in5(:,j));
    OA(1,j) = 100*sum(c==prd)/length(c);
end
[a1,b1]=find(OA==max(OA(:)));
a1=a1(1);b1=b1(1);
disp(max(OA(:)))
disp(a(a1))
disp(b1)
gmm=fitgmdist(in5(:,b1),2);
P = posterior(gmm,(-1:.01:1)');
[~,w]=(min(abs(P(:,1)-0.5)));
ThresholdG0=-1+w*0.01;
[~,~,c]=unique(cl5);
c(c==a1)=10;
c(c~=10)=11;
mld = fitctree(in5(:,b1),c,'MaxNumSplits',1);
view(mld,'Mode','graph')
prd = mld.predict(in5(:,b1));
OA= 100*sum(c==prd)/length(c);

%% Training results based on above

% Windices=[NDVI,MNDWI1,NDWI,NDTI,NDMI1,...
%     EVI,FAI,AFAI,ABDI,TMI,...
%     SABI,VBFAH,SA,MAI,SUI,...
%     NDCI,MOSES3B,S23BDA,B3B2,GNDVI,...
%     EVSI,CMI,NDAVI,WAVI];


clearvars -except cl25 in25
cl5=cl25;
in5=in25;

[~,~,c]=unique(cl5);
pr=nan(size(cl5));
idx = find(in5(:,14)>0.152); % MAI
pr(idx) = 3; % Land
in5(idx,:)=nan;
idx = find(in5(:,8)<-0.033); % AFAI
pr(idx) = 1; % Algae
in5(idx,:)=nan;
idx = find(in5(:,13)>-0.1); % SA
pr(idx) = 5; % Water
in5(idx,:)=nan;
idx = in5(:,7)>0.046; % FAI
pr(idx) = 2; % Emergent
pr(isnan(pr)) = 4; % Submerged

OA = 100*sum(c==pr)/length(c);

figure;plotconfusion(categorical(c),categorical(pr))


figure;
subplot(2,4,1);boxplot(in25(:,17),cl25(:));title('MAI');
subplot(2,4,2);boxplot(in25(:,10),cl25(:));title('AFAI');
subplot(2,4,3);boxplot(in25(:,9),cl25(:));title('FAI');
subplot(2,4,4);boxplot(in25(:,16),cl25(:));title('SA');
subplot(2,4,4+1);boxplot(in25(:,3),cl25(:));title('MNDWI2');
subplot(2,4,4+2);boxplot(in25(:,11),cl25(:));title('ABDI');
subplot(2,4,4+3);boxplot(in25(:,14),cl25(:));title('IGAG');
subplot(2,4,4+4);boxplot(in25(:,22),cl25(:));title('Toming');

figure
subplot(1,4,1);boxplot(in25(:,3),cl25(:),"Colors",'bkgcr','Widths',0.8,Whisker=0,OutlierSize=0.02);title('MNDWI2');
subplot(1,4,2);boxplot(in25(:,1),cl25(:),"Colors",'bkgcr','Widths',0.8,Whisker=0,OutlierSize=0.02);title('NDVI');
subplot(1,4,3);boxplot(in25(:,26),cl25(:),"Colors",'bkgcr','Widths',0.8,Whisker=0,OutlierSize=0.02);title('GNDVI');
subplot(1,4,4);boxplot(in25(:,17),cl25(:),"Colors",'bkgcr','Widths',0.8,Whisker=0,OutlierSize=0.02);title('MAI');

figure
subplot(1,4,1);boxplot(in25(:,17),cl25(:),"Colors",'bkgcr','Widths',0.8,Whisker=0,OutlierSize=0.02);title('MAI');
subplot(1,4,2);boxplot(in25(:,10),cl25(:),"Colors",'bkgcr','Widths',0.8,Whisker=0,OutlierSize=0.02);title('AFAI');
subplot(1,4,3);boxplot(in25(:,16),cl25(:),"Colors",'bkgcr','Widths',0.8,Whisker=0,OutlierSize=0.02);title('SA');
subplot(1,4,4);boxplot(in25(:,9),cl25(:),"Colors",'bkgcr','Widths',0.8,Whisker=0,OutlierSize=0.02);title('FAI');



clearvars -except cl25 in25
cl5=cl25;
in5=in25;
[~,~,c]=unique(cl5);
pr=nan(size(cl5));
idx = find(in5(:,17)>0.150); % MAI
pr(idx) = 3; % Land
in5(idx,:)=nan;
idx = find(in5(:,10)<-0.03); % AFAI
pr(idx) = 1; % Algae
in5(idx,:)=nan;
idx = find(in5(:,16)>-0.05); % SA
pr(idx) = 5; % Water
in5(idx,:)=nan;
idx = in5(:,9)>0.045; % FAI
pr(idx) = 2; % Emergent
pr(isnan(pr)) = 4; % Submerged
OA = 100*sum(c==pr)/length(c);
figure;plotconfusion(categorical(c),categorical(pr))

c(c==4)=2;
pr(pr==4)=2;
figure;plotconfusion(categorical(c),categorical(pr))









in04=[in25(:,9),in25(:,10),in25(:,16),in25(:,17)];


% prTr=mld.predict(in04);
% figure;plotconfusion(categorical(cl5),categorical(prTr))


%% Test on new dataset
clearvars -except cl25 in25

data1=readcell("G:\Trainee\pragati\Nal_100_DataPoints.xlsx");
ref1=cell2mat(data1(2:end,5:16));
class1=data1(2:end,25);
ref=ref1./10000;
Wclass=class1;

FAI =  ref(:,8) - (ref(:,4) + (ref(:,11)-ref(:,4))*0.1873); % Floating Algae Index
AFAI = ref(:,8) - ref(:,4) + (ref(:,11)-ref(:,4))*0.5; % Alternative Floating Algae Index
SA = (ref(:,4)-ref(:,8))/.177; % Slope Algorithm
MAI = ref(:,3) + ref(:,4) - (ref(:,11)+(ref(:,2)-ref(:,11))*0.9375); % Macro-Algae Index

in4=[FAI,AFAI,SA,MAI];

[~,~,c]=unique(Wclass);
pr=nan(size(Wclass));
idx = find(in4(:,4)>0.152);
pr(idx) = 3; % Land
in4(idx,:)=nan;
idx = find(in4(:,2)<-0.033);
pr(idx) = 1; % Algae
in4(idx,:)=nan;
idx = find(in4(:,3)>-0.1);
pr(idx) = 5;  % Water
in4(idx,:)=nan;
idx = in4(:,1)>0.046;
pr(idx) = 2; % Emergent
pr(isnan(pr)) = 4;  % Submerged

OA = 100*sum(c==pr)/length(c);

figure;plotconfusion(categorical(c),categorical(pr))


c(c==4)=2;
pr(pr==4)=2;
figure;plotconfusion(categorical(c),categorical(pr))





mld = fitctree(in4,Wclass,"MaxNumSplits",5);
view(mld,'Mode','graph')
prTr=mld.predict(in4);
figure;plotconfusion(categorical(Wclass),categorical(prTr))



figure;
subplot(1,4,1);boxplot(in4(:,4),Wclass);title('MAI');
subplot(1,4,2);boxplot(in4(:,2),Wclass);title('AFAI');
subplot(1,4,3);boxplot(in4(:,3),Wclass);title('SA');
subplot(1,4,4);boxplot(in4(:,1),Wclass);title('FAI');

%% Test on krishna dataset
clearvars -except cl25 in25 a c

data1=readcell("G:\Trainee\pragati\Sample_Classification_krishna.xlsx");
ref1=cell2mat(data1(2:end,5:16));
% ref1([26,35,36],:)=[];
class1=data1(2:end,25);
% class1([26,35,36])=[];
ref=ref1./10000;
Wclass=class1;

FAI =  ref(:,8) - (ref(:,4) + (ref(:,11)-ref(:,4))*0.1873); % Floating Algae Index
AFAI = ref(:,8) - ref(:,4) + (ref(:,11)-ref(:,4))*0.5; % Alternative Floating Algae Index
SA = (ref(:,4)-ref(:,8))/.177; % Slope Algorithm
MAI = ref(:,3) + ref(:,4) - (ref(:,11)+(ref(:,2)-ref(:,11))*0.9375); % Macro-Algae Index

in4=[FAI,AFAI,SA,MAI];

[~,~,c]=unique(Wclass);
pr=nan(size(Wclass));
idx = find(in4(:,4)>0.152);
pr(idx) = 3; % Land
in4(idx,:)=nan;
idx = find(in4(:,2)<-0.033);
pr(idx) = 1; % Algae
in4(idx,:)=nan;
idx = find(in4(:,3)>-0.1);
pr(idx) = 5;  % Water
in4(idx,:)=nan;
idx = in4(:,1)>0.046;
pr(idx) = 2; % Emergent
pr(isnan(pr)) = 4;  % Submerged

OA = 100*sum(c==pr)/length(c);

figure;plotconfusion(categorical(c),categorical(pr))


c(c==4)=1;
pr(pr==4)=1;
figure;plotconfusion(categorical(c),categorical(pr))

figure;
subplot(1,4,1);boxplot(in4(:,4),Wclass);title('MAI');
subplot(1,4,2);boxplot(in4(:,2),Wclass);title('AFAI');
subplot(1,4,3);boxplot(in4(:,3),Wclass);title('SA');
subplot(1,4,4);boxplot(in4(:,1),Wclass);title('FAI');



NDVI=(ref(:,8)-ref(:,4))./(ref(:,8)+ref(:,4)); % Normalized Difference Vegetation Index
MNDWI2=(ref(:,3)-ref(:,12))./(ref(:,3)+ref(:,12)); % Modified Normalized Difference Water Index using SWIR band 2
GNDVI = (ref(:,8)-ref(:,3))./(ref(:,8)+ref(:,3)); % Green Normalized Difference Vegetation Index
AFAI = ref(:,8) - ref(:,4) + (ref(:,11)-ref(:,4))*0.5; % Alternative Floating Algae Index

in4=[MNDWI2,NDVI,AFAI,GNDVI];


[~,~,c]=unique(Wclass);
pr=nan(size(Wclass));
idx = find(in4(:,1)<-0.15); % MNDWI2
pr(idx) = 3; % Land
in4(idx,:)=nan;
idx = find(in4(:,2)<0.1 & in4(:,3)>-0.025); % NDVI and AFAI
pr(idx) = 5; % water
in4(idx,:)=nan;
idx = find(in4(:,2)>0.35); % NDVI
pr(idx) = 2; % emergent
in4(idx,:)=nan;
idx = in4(:,4)<0; % GNDVI
pr(idx) = 1; % Algae
pr(isnan(pr)) = 4; % Submerged
OA = 100*sum(c==pr)/length(c);
figure;plotconfusion(categorical(c),categorical(pr))


figure;
subplot(1,4,1);boxplot(in4(:,4),Wclass);title('GNDVI');
subplot(1,4,2);boxplot(in4(:,2),Wclass);title('NDVI');
subplot(1,4,3);boxplot(in4(:,3),Wclass);title('AFAI');
subplot(1,4,4);boxplot(in4(:,1),Wclass);title('MNDWI2');





%%
% clearvars -except cl25 in25
% 
% data1=readcell("G:\Trainee\pragati\Ukai and Vembanad data points.xlsx");
% ref1=cell2mat(data1(2:end,5:16));
% class1=data1(2:end,25);
% ref=ref1./10000;
% Wclass=class1;
% 
% FAI =  ref(:,8) - (ref(:,4) + (ref(:,11)-ref(:,4))*0.1873); % Floating Algae Index
% AFAI = ref(:,8) - ref(:,4) + (ref(:,11)-ref(:,4))*0.5; % Alternative Floating Algae Index
% SA = (ref(:,4)-ref(:,8))/.177; % Slope Algorithm
% MAI = ref(:,3) + ref(:,4) - (ref(:,11)+(ref(:,2)-ref(:,11))*0.9375); % Macro-Algae Index
% 
% in4=[FAI,AFAI,SA,MAI];
% 
% [~,~,c]=unique(Wclass);
% pr=nan(size(Wclass));
% idx = find(in4(:,4)>0.152); % MAI
% pr(idx) = 3; % Land
% in4(idx,:)=nan;
% idx = find(in4(:,2)<-0.033); % AFAI
% pr(idx) = 1; % Algae
% in4(idx,:)=nan;
% idx = find(in4(:,3)>-0.1); % SA
% pr(idx) = 5;  % Water
% in4(idx,:)=nan;
% idx = in4(:,1)>0.046; % FAI
% pr(idx) = 2; % Emergent
% pr(isnan(pr)) = 4;  % Submerged
% 
% OA = 100*sum(c==pr)/length(c);
% 
% figure;plotconfusion(categorical(c),categorical(pr))

%%


% clearvars -except cl25 in25
% cl5=cl25;
% in5=in25;
% [a,~,c]=unique(cl5);
% pr=nan(size(cl5));
% idx = find(in5(:,17)>0.150 & in5(:,3)<-0.15); % MAI & MNDWI2
% pr(idx) = 3; % Land
% in5(idx,:)=nan;
% idx = find(in5(:,10)<-0.03 & in5(:,11)<-0.03); % AFAI & ABDI
% pr(idx) = 1; % Algae
% in5(idx,:)=nan;
% idx = find(in5(:,9)<-0.01 & in5(:,16)>-0.1); % FAI
% pr(idx) = 5; % Water
% in5(idx,:)=nan;
% idx = find(in5(:,16)<-0.3 & in5(:,22)<0.07); % SA
% pr(idx) = 2; % Emergent
% pr(isnan(pr)) = 4; % Submerged
% OA = 100*sum(c==pr)/length(c);
% figure;plotconfusion(categorical(c),categorical(pr))


clearvars -except cl25 in25
cl5=cl25;
in5=in25;
[~,~,c]=unique(cl5);
pr=nan(size(cl5));
idx = find(in5(:,17)>0.150); % MAI
pr(idx) = 3; % Land
in5(idx,:)=nan;
idx = find(in5(:,10)<-0.03); % AFAI
pr(idx) = 1; % Algae
in5(idx,:)=nan;
idx = find(in5(:,16)>-0.05); % SA
pr(idx) = 5; % Water
in5(idx,:)=nan;
idx = in5(:,9)>0.045; % FAI
pr(idx) = 2; % Emergent
pr(isnan(pr)) = 4; % Submerged
OA = 100*sum(c==pr)/length(c);
figure;plotconfusion(categorical(c),categorical(pr))


clearvars -except cl25 in25
cl5=cl25;
in5=in25;
[a,~,c]=unique(cl5);
pr=nan(size(cl5));
idx = find(in5(:,3)<-0.15); % MNDWI2
pr(idx) = 3; % Land
in5(idx,:)=nan;
idx = find(in5(:,1)<0.1 & in5(:,10)>-0.025); % NDVI and AFAI
pr(idx) = 5; % water
in5(idx,:)=nan;
idx = find(in5(:,1)>0.35); % NDVI
pr(idx) = 2; % emergent
in5(idx,:)=nan;
idx = find(in5(:,26)<0); % GNDVI
pr(idx) = 1; % Algae
pr(isnan(pr)) = 4; % Submerged
OA = 100*sum(c==pr)/length(c);
figure;plotconfusion(categorical(c),categorical(pr))


[~,~,c]=unique(cl25);
figure;
x=-1:0.001:1;
for i=1:30
    if i==14
    else
    subplot(5,6,i);
    hold on
    for j=1:5
    mld=fitdist(in25(c==j,i),"Kernel");
    y=mld.pdf(x);
    p=area(x,y);
    p(1).FaceAlpha=0.5;
    end
    hold off
    end
end


