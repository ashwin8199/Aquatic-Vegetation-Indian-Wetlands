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

in25(in25>100)=0;
in25(in25<-100)=0;
clearvars -except cl25 in25 



data1=readcell("G:\Trainee\pragati\Sample_Classification_krishna.xlsx");
ref1=cell2mat(data1(2:end,5:16));
% ref1([26,35,36],:)=[];
class1=data1(2:end,25);
% class1([26,35,36])=[];
ref=ref1./10000;
Wclass=class1;

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


Tcl25=Wclass(maj5in);
Tin25= Windices(maj5in,:);
Tin25(Tin25>100)=0;
Tin25(Tin25<-100)=0;

clearvars -except cl25 in25 maj5in Tcl25 Tin25

%%



rfmodel = TreeBagger(50,in25,cl25,NumPredictorsToSample=24,OOBPrediction="on",Method="classification");
oobErr= oobError(rfmodel);
oobAcc = 1-oobErr(end)
Ypred = predict(rfmodel,Tin25);
testACC = sum(strcmp(Ypred,Tcl25))/length(Tcl25)
view(rfmodel.Trees{1},'Mode','graph')

