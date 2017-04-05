clear all
close
clc

%Location of files
[~, name] = system('hostname'); 

folder='K:\windows\data\RoyData\160205_Agar_Data';
if strcmp( name, 'Atlantis') %Josko home PC
    folder = 'K:\windows\data\RoyData\160205_Agar_Data';
end
slash = '/';
%read experiments sequentially
exps=[1 2 3 4 5 7 8 9];
Intensityval = [2000, 1500, 3300].*2; %[CFP YFP RFP]
umperpx=0.159;
%For each channel
%   Length, Long-Axis Position, Integrated Intensity,, Full-Cell Intensity,
%       number of spots, cell lengths
Lcfp=[];    Lyfp=[];    Lrfp=[];
Pcfp=[];    Pyfp=[];    Prfp=[];
Icfp=[];    Iyfp=[];    Irfp=[];
Fcfp=[];    Fyfp=[];    Frfp=[];
ncfp=[];    nyfp=[];    nrfp=[];
celllength = [];

%   plotting parameters
yfp.filterval=Intensityval(2)*2;
cfp.filterval=Intensityval(1)*2;
rfp.filterval=Intensityval(3)*2;

%   pure loading
j=1; 
for i=exps;
    E{j}=load(strcat(folder,slash,num2str(i),slash,'Results.mat')); 
    j=j+1;
end

allCFP_L = [];

Nexp=size(E,2);

%reformat data
for i=1:Nexp

    Ncells{i}=size(E{i}.DataStruct,2); 
    for j=1:Ncells{i} 

        if ~isempty(E{i}.DataStruct(1,j).Lnorm)        
            LNormCFP{i,j}=E{i}.DataStruct(1,j).Lnorm;
        else
            LNormCFP{i,j}=0;
        end

        CellLength{i,j}=E{i}.DataStruct(1,j).CellLength;


        CFPld{i,j}=E{i}.DataStruct(1,j).ld;
        YFPld{i,j}=E{i}.DataStruct(2,j).ld;
        RFPld{i,j}=E{i}.DataStruct(3,j).ld;

        NspotsCFP=size(CFPld{i,j},2);
        NspotsYFP=size(YFPld{i,j},2);
        NspotsRFP=size(RFPld{i,j},2);        

        if NspotsCFP==0
            CFPld{i,j}{1}=[];
        else
            for k=1:NspotsCFP
                Lcfp=[Lcfp CellLength{i,j}];
                Pcfp=[Pcfp CFPld{i,j}{k}(1,2)/CellLength{i,j}];
                Icfp=[Icfp 2*pi*CFPld{i,j}{k}(1,1)*CFPld{i,j}{k}(1,3)*CFPld{i,j}{k}(1,5)/Intensityval(1)];
                Fcfp=[Fcfp CFPld{i,j}{k}(1,7)];
            end
        end
        ncfp = [ncfp NspotsCFP];
        celllength = [celllength CellLength{i,j}];

        if NspotsYFP==0
            YFPld{i,j}{1}=[];
        else
            for k=1:NspotsYFP
                Lyfp=[Lyfp CellLength{i,j}];
                Pyfp=[Pyfp YFPld{i,j}{k}(1,2)/CellLength{i,j}];
                Iyfp=[Iyfp 2*pi*YFPld{i,j}{k}(1,1)*YFPld{i,j}{k}(1,3)*YFPld{i,j}{k}(1,5)/Intensityval(2)];
                Fyfp=[Fyfp YFPld{i,j}{k}(1,7)];
            end
        end
        nyfp = [nyfp NspotsYFP];

        if NspotsRFP==0
            RFPld{i,j}{1}=[];
        else
            for k=1:NspotsRFP            
                Lrfp=[Lrfp CellLength{i,j}];
                Prfp=[Prfp RFPld{i,j}{k}(1,2)/CellLength{i,j}];
                Irfp=[Irfp 2*pi*RFPld{i,j}{k}(1,1)*RFPld{i,j}{k}(1,3)*RFPld{i,j}{k}(1,5)/Intensityval(3)];
                Frfp=[Frfp RFPld{i,j}{k}(1,7)];
            end
        end
        nrfp = [nrfp NspotsRFP];
    end
end

fprintf('Data has been loaded and re-formatted. Next, figures.\n'); 
 

%Simple figure showing scattered spot positions versus cell length. Note that cell length is
%   related to but not equal to the replication cycle, causing a large error both vertically
%   and horizontally.

%   The markers are at constant opacity, to give an idea of density, and of a marker size that is
%   proportional to their relative intensity. 

opacity = 0.6;
 
fig1 = figure(1);  
set(fig1,'Position',[100,0,1820,1080])
subplot(1,3,1);
spotPositionCellLength(Lcfp, Pcfp, Icfp, 'CFP', 'b', opacity);
set(gca,'Color',[0. 0. 0.]);

subplot(1,3,2);
spotPositionCellLength(Lyfp, Pyfp, Iyfp, 'YFP', 'y', opacity);
set(gca,'Color',[0. 0. 0.]);

subplot(1,3,3);
spotPositionCellLength(Lrfp, Prfp, Irfp, 'RFP', 'r', opacity);
set(gca,'Color',[0. 0. 0.]);

%   This figure shows the spot position versus intensity. THe goal of this is to see if teh
%cell length is correlated to the spot intensity, but the results so far clearly show that this 
%is mostly random. A polynomial fit doesn't tell us much; you need a model before doing so.


fig2 = figure(2);
set(fig2,'Position',[100,0,1820,1080])
subplot(1,3,1) 
spotPositionIntensity(Pcfp, Icfp, 'b', 'CFP');
subplot(1,3,2) 
spotPositionIntensity(Pyfp, Iyfp, 'y', 'YFP');
subplot(1,3,3) 
spotPositionIntensity(Prfp, Irfp, 'r', 'RFP');
 
 
% This shows the spots again, but with histograms showing the distribution of spots.
%       bin size is calculated by Sturges' formula
fig3 = figure(3);
set(fig3,'Position',[100,0,1820,1080]) 

subplot(1,3,1);
spotPositionCount(Pcfp, Icfp, 'b', 'CFP');
set(gca,'Color',[0. 0. 0.]);
subplot(1,3,2);
spotPositionCount(Pyfp, Iyfp, 'y', 'YFP');
set(gca,'Color',[0. 0. 0.]);
subplot(1,3,3);
spotPositionCount(Prfp, Irfp, 'r', 'RFP'); 
set(gca,'Color',[0. 0. 0.]);

abort('Figure %d', 3);
bins = 20;
thisedge2{1} = linspace(15,35,bins+1)*0.159;
thisedge2{2} = (0:bins)/bins;
 
set(fig4,'Position',[20,300,1800,500]) 
set(fig5,'Position',[20,300,1800,500])

% CFP

subplot(1,3,1)
Numcfp(1,:) = Lcfp*0.159;
Numcfp(2,:) = Pcfp;

%Filter on Intensity
ICFP=Icfp*Intensityval(1);
[bin_cfp,idx_cfp]=find(ICFP>cfp.filterval); 

NumCFP(:,idx_cfp)=Numcfp(:,idx_cfp);
NumCFPnz(1,:)=nonzeros(NumCFP(1,:));
NumCFPnz(2,:)=nonzeros(NumCFP(2,:));

FilteredSignalCFP=size(NumCFPnz,2)/size(ICFP,2); % Percentage of total signal

pause;figure(4)
subplot(1,3,1)
Heatmap = hist3(NumCFPnz','Edges',thisedge2);
Heatmap=Heatmap';
DummyHeat=max(Heatmap);
for i=1:size(Heatmap,2)
Heatmap(:,i)=(Heatmap(:,i)./DummyHeat(i));
end
pcolor(thisedge2{1},(thisedge2{2}),Heatmap);
% colormap(fig4,jet) % heat map
xlabel('Cell Length'); ylabel('Position in Cell');
title('Agar data: CFP');
grid on
set(gca,'FontSize',16)

pause;figure(5)
subplot(1,3,1)
hold on
hist3(NumCFPnz','Edges',thisedge2)
colormap(fig5,jet) % heat map
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
xlabel('Cell Length'); ylabel('Position in Cell');zlabel('Amount of spots')
title('Agar data: CFP');
grid off
hold off
axis([min(Lcfp),max(Lcfp),0,1])
view(3)
set(gca,'FontSize',16)

% YFP
Numyfp(1,:) = Lyfp*0.159;
Numyfp(2,:) = Pyfp;

%Filter on Intensity
IYFP=Iyfp*Intensityval(2);
[bin_yfp,idx_yfp]=find(IYFP>yfp.filterval); 

NumYFP(:,idx_yfp)=Numyfp(:,idx_yfp);
NumYFPnz(1,:)=nonzeros(NumYFP(1,:));
NumYFPnz(2,:)=nonzeros(NumYFP(2,:));

FilteredSignalYFP=size(NumYFPnz,2)/size(IYFP,2); % Percentage of total signal

pause;figure(4)
subplot(1,3,2)
Heatmap = hist3(NumYFPnz','Edges',thisedge2);
Heatmap=Heatmap';
DummyHeat=max(Heatmap);
for i=1:size(Heatmap,2)
Heatmap(:,i)=(Heatmap(:,i)./DummyHeat(i));
end
pcolor(thisedge2{1},(thisedge2{2}),Heatmap);
xlabel('Cell Length'); ylabel('Position in Cell');
title('Agar data: YFP');
grid on
set(gca,'FontSize',16)

pause;figure(5)
subplot(1,3,2)
hold on
hist3(NumYFPnz','Edges',thisedge2)
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
xlabel('Cell Length'); ylabel('Position in Cell');zlabel('Amount of spots')
title('Agar data: YFP');
grid off
hold off
axis([min(Lyfp),max(Lyfp),0,1])
view(3)
set(gca,'FontSize',16)

% RFP
Numrfp(1,:) = Lrfp*0.159;
Numrfp(2,:) = Prfp;

%Filter on Intensity
IRFP=Irfp*Intensityval(3);
[bin_rfp,idx_rfp]=find(IRFP>rfp.filterval); 

NumRFP(:,idx_rfp)=Numrfp(:,idx_rfp);

NumRFPnz(1,:)=nonzeros(NumRFP(1,:));
NumRFPnz(2,:)=nonzeros(NumRFP(2,:));

FilteredSignalRFP=size(NumRFPnz,2)/size(IRFP,2);

pause;figure(4)
subplot(1,3,3)
Heatmap = hist3(NumRFPnz','Edges',thisedge2);
Heatmap=Heatmap';
DummyHeat=max(Heatmap);
for i=1:size(Heatmap,2)
Heatmap(:,i)=(Heatmap(:,i)./DummyHeat(i));
end
h = pcolor(thisedge2{1},(thisedge2{2}),Heatmap);
xlabel('Cell Length'); ylabel('Position in Cell');
title('Agar data: RFP');
grid on
set(gca,'FontSize',16)

pause;figure(5)
subplot(1,3,3)
hold on
hist3(NumRFPnz','Edges',thisedge2)
colormap(jet) % heat map
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
xlabel('Cell Length'); ylabel('Position in Cell');zlabel('Amount of spots')
title('Agar data: RFP');
grid off
hold off
axis([min(Lrfp),max(Lrfp),0,1])
view(3)
set(gca,'FontSize',16)

%% Full cell intensity vs. celllength

clear plotcfp plotyfp plotrfp
pause;fig6 = figure(6);
set(fig6,'Position',[20,300,1800,500])

% CFP
bins1 = 7;
bins2 = 50;
thisedge3{1} = linspace(min(Lcfp),max(Lcfp),bins1+1);
thisedge3{2} = linspace(min(Fcfp),max(Fcfp),bins2+1);

plotcfp(1,:) = Lcfp; %length
plotcfp(2,:) = Fcfp; %full cell intensity
plotcfp(3,:) = Icfp*Intensityval(1); %spot intensity

plotcfp = unique(plotcfp','rows')';
% 
% [N_full,Edges_full,mid_full,loc_full]=histcn([plotcfp(1,:)' plotcfp(2,:)'],thisedge3{1},thisedge3{2});
% [N_spot,Edges_spot,mid_spot,loc_spot]=histcn([plotcfp(1,:)' plotcfp(3,:)'],thisedge3{1},thisedge3{2});
% 
% for i=1:size(N_full,1);
% CFP_binned(i,:)=N_full(i,:).*Edges_full{2};
% Y=CFP_binned(i,:)';
% Y=double(Y);
% X=1:length(CFP_binned(i,:));
% display(strcat({'Gaussian fit of CFP column '},num2str(i),{' of '},num2str(size(N_full,1))));
% f{i}=fit(X',Y,'gauss1');
% peak(i)=f{i}.b1;
% sigma(i)=f{i}.c1;
% peakfloored(i)=floor(peak(i));
% peakceiled(i)=ceil(peak(i));
% IntVal(i)=((Edges_full{2}(peakceiled(i))-Edges_full{2}(peakfloored(i))))*(peakceiled(i)-peak(i))+Edges_full{2}(peakfloored(i)); %slope * peakposition, because linear.
% end


subplot(1,3,1)
hold on
scatter(plotcfp(1,:),plotcfp(2,:),'b','o','filled');
scatter(plotcfp(1,:),plotcfp(3,:),'r','o','filled');
myfit=polyfit(plotcfp(1,:),plotcfp(2,:),1);
myfit2=polyfit(plotcfp(1,:),plotcfp(3,:),1);
x=12:0.1:43;
y=polyval(myfit,x);
y2=polyval(myfit2,x);
plot(x,y,'b','LineWidth',3)
plot(x,y2,'r','LineWidth',3)
xlabel('Cell Length (px)'); ylabel('Intensity (-)'); 
title('TetR stoichiometry vs. Length')
hold off
axis([12 43 -0.1 2*10^5])
set(gca,'FontSize',16)
legend('Full cell intensity','Spot intensity')

clear N_full
% YFP


plotyfp(1,:) = Lyfp;
plotyfp(2,:) = Fyfp;
plotyfp(3,:) = Iyfp;

plotyfp = unique(plotyfp','rows')';

[N_full,Edges_full,mid_full,loc_full]=histcn([plotyfp(1,:)' plotyfp(2,:)'],thisedge3{1},thisedge3{2});
[N_spot,Edges_spot,mid_spot,loc_spot]=histcn([plotyfp(1,:)' plotyfp(3,:)'],thisedge3{1},thisedge3{2});

% for i=1:size(N_full,1);
% YFP_binned(i,:)=N_full(i,:);
% Y=YFP_binned(i,:)';
% Y=double(Y);
% X=1:length(YFP_binned(i,:));
% display(strcat({'Gaussian fit of CFP column '},num2str(i),{' of '},num2str(size(N_full,1))));
% f{i}=fit(X',Y,'gauss1');
% peak(i)=f{i}.b1;
% sigma(i)=f{i}.c1;
% peakfloored(i)=floor(peak(i));
% peakceiled(i)=ceil(peak(i));
% IntVal(i)=((Edges_full{2}(peakceiled(i))-Edges_full{2}(peakfloored(i))))*(peakceiled(i)-peak(i))+Edges_full{2}(peakfloored(i)); %slope * peakposition, because linear.
% end
% plot(Edges_full{1},IntVal);

subplot(1,3,2)
hold on
scatter(plotyfp(1,:),plotyfp(2,:)/Intensityval(2),'b','o','filled');
scatter(plotyfp(1,:),plotyfp(3,:)/Intensityval(2),'r','o','filled');
myfit=polyfit(plotyfp(1,:),plotyfp(2,:)/Intensityval(2),1);
myfit2=polyfit(plotyfp(1,:),plotyfp(3,:)/Intensityval(2),1);
x=12:0.1:43;
y=polyval(myfit,x);
y2=polyval(myfit2,x);
plot(x,y,'b','LineWidth',3)
% plot(x,y2,'r','LineWidth',3)
xlabel('Cell Length'); ylabel('Intensity'); 
title('Tus Stoichiometry vs. Length')
hold off
% axis([12 43 -0.1 1])
set(gca,'FontSize',16)


% RFP

plotrfp(1,:) = Lrfp;
plotrfp(2,:) = Frfp;
plotrfp(3,:) = Irfp*Intensityval(3);
plotrfp = unique(plotrfp','rows')';

subplot(1,3,3)
hold on
scatter(plotrfp(1,:),plotrfp(2,:),'b','o','filled');
scatter(plotrfp(1,:),plotrfp(3,:),'r','o','filled');
myfit=polyfit(plotrfp(1,:),plotrfp(2,:),1);
myfit2=polyfit(plotrfp(1,:),plotrfp(3,:),1);
x=12:0.1:43;
y=polyval(myfit,x);
y2=polyval(myfit2,x);
plot(x,y,'b','LineWidth',3)
plot(x,y2,'r','LineWidth',3)
xlabel('Cell Length'); ylabel('Normalized full cell intensity'); 
title('DnaN Stoichiometry vs. Length')
hold off
axis([12 43 -0.1 4*10^5])
set(gca,'FontSize',16)

%% Numspots/cell vs. cell length
pause;fig1 = figure(7);

Datapoints = 6; %number of points for the average


for i=1:max(nrfp)%i=1:Datapoints;
p{i} = find(nrfp == i);
m(i) = mean(celllength(p{i}));
Spotstd(i)=std(celllength(p{i}));
%Interval=((max(celllength)-min(celllength))/Datapoints)+i*min(celllength);
%p{i}=find(celllength >= 
end

myfit=polyfit(m,1:4,1);
x=[12,43];
y=polyval(myfit,x);

Y=linspace(1,max(nrfp),max(nrfp));

Xl=m-Spotstd;
Xr=m+Spotstd;

hold on
scatter(celllength,nrfp,'m')
plot(m,Y,'b--','LineWidth',3)
plot(Xl,Y,'b','LineWidth',1)
plot(Xr,Y,'b','LineWidth',1)
axis([12, 43, 0, 5])
xlabel('Cell length'); ylabel('Number of spots per cell')
title('Numspots/cell vs. cell length for YFP')
hold off
set(gca,'FontSize',16)
legend('Spot vs Length points','Mean value','std')

pause;
close all;