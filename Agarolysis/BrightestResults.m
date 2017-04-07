clear all
close
clc

%Location of files
[~, name] = system('hostname'); 

folder='D:\RoyData\160205_Agar_Data';
if strcmp( name, 'Atlantis') %Josko home PC
    folder = 'K:\windows\data\RoyData\160205_Agar_Data';
end
slash = '/';
%read experiments sequentially
exps=[1 2 3 4 5 7 8 9];
calibration = [2000, 1500, 3300].*2; %[CFP YFP RFP]
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

%   pure loading
j=1; 
for i=exps;
    E{j}=load(strcat(folder,slash,num2str(i),slash,'Results.mat')); 
    j=j+1;
end

allCFP_L = [];

Nexp=size(E,2);

pairwiseLength = [];
pairwiseC = [];
pairwiseY = [];
pairwiseR = [];

%reformat data
for i=1:Nexp

    Ncells{i}=size(E{i}.DataStruct,2); 
    for j=1:Ncells{i} 
        fprintf('Experiment %d, cell %d. \n', i, j);
        if ~isempty(E{i}.DataStruct(1,j).Lnorm)        
            LNormCFP{i,j}=E{i}.DataStruct(1,j).Lnorm;
        else
            LNormCFP{i,j}=0;
        end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        skipPair = 0;
        if length(E{i}.DataStruct(1,j).ld) == 0
            skipPair = 1;
        end
        if length(E{i}.DataStruct(2,j).ld) == 0
            skipPair = 1;
        end
        if length(E{i}.DataStruct(3,j).ld) == 0 
            skipPair = 1;
        end
        if skipPair == 0

            fig1 = figure(1);
            set(fig1,'Position',[0,0,1680,1080])
            for k=1:3
                channel = 'CFP';
                if k == 2
                    channel = 'YFP';
                end
                if k == 3
                    channel = 'RFP' ;
                end
                fprintf('Experiment %d, Cell %d, channel %s: \n', i,j, k);

                subplot(1,3,k);  
                hold on;
                    surf( E{i}.DataStruct(k,j).ydatacrpdR1{1}); 

                hold off;
                colormap hot;
                shading interp;
                view([0 90])

                fitData = E{i}.DataStruct(k,j).x{:};
                fprintf('\tSpot %.3f, %.3f+-%.3f, %.3f+-%.3f \n', fitData(1), fitData(2), fitData(3), fitData(4), fitData(5)) 

                badSpots = E{i}.DataStruct(3,j).bx;
                fprintf('\tRejected:\n');

                for q = 1:length(badSpots)
                    spot = badSpots{q}; 
                    fprintf('\tSpot %.3f, %.3f+-%.3f, %.3f+-%.3f \n', spot(1), spot(2), spot(3), spot(4), spot(5)) 
                end 
            end

            fitData = E{i}.DataStruct(1,j).ld{:}; 
            pairwiseC = [pairwiseC, [fitData(2) fitData(4)]'];

            fitData = E{i}.DataStruct(2,j).ld{:};
            pairwiseY = [pairwiseY, [fitData(2) fitData(4)]'];

            fitData = E{i}.DataStruct(3,j).ld{:};
            pairwiseR = [pairwiseR, [fitData(2) fitData(4)]'];

            pairwiseLength = [pairwiseLength, E{i}.DataStruct(1,j).CellLength]; 

            pause; 
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
                Icfp=[Icfp 2*pi*CFPld{i,j}{k}(1,1)*CFPld{i,j}{k}(1,3)*CFPld{i,j}{k}(1,5)/calibration(1)];
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
                Iyfp=[Iyfp 2*pi*YFPld{i,j}{k}(1,1)*YFPld{i,j}{k}(1,3)*YFPld{i,j}{k}(1,5)/calibration(2)];
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
                Irfp=[Irfp 2*pi*RFPld{i,j}{k}(1,1)*RFPld{i,j}{k}(1,3)*RFPld{i,j}{k}(1,5)/calibration(3)];
                Frfp=[Frfp RFPld{i,j}{k}(1,7)];
            end
        end
        nrfp = [nrfp NspotsRFP];
    end
end
close all;
fprintf('Data has been loaded and re-formatted. Next, figures.\n'); 
figNum = 0;
if 0
    %scatter [cell length, distance between brightest points for |YFP-CFP| (Tus-Dif), |YFP-RFP| (TUS-DNAN) and |YFP-RAND| (TUS random)]

    lengthCR = sqrt( sum((pairwiseC - pairwiseR).^2));
    lengthCY = sqrt( sum((pairwiseC - pairwiseY).^2));
    lengthRY = sqrt( sum((pairwiseR - pairwiseY).^2));

    scatter(pairwiseLength, lengthCR);
    pause;
    scatter(pairwiseLength, lengthCY);
    pause;
    scatter(pairwiseLength, lengthRY);


    pause;close all; clc;
    error('abort');
end
%error(' abort ' );
figNum = figNum + 1; 
fig0 = figure(figNum);
histogram( Lcfp, linspace( floor(min(Lcfp)), ceil(max(Lcfp)), log(length(Lcfp))/log(2)+1));
xlabel('Cell Length');
ylabel('Counts');
title('Select left limit for scatter');


[leftLim, dummy] = ginput(1);
title('Select right limit for scatter');
[rightLim, dummy] = ginput(1);

%abort('Hey %s', 'test');

%Simple figure showing scattered spot positions versus cell length. Note that cell length is
%   related to but not equal to the replication cycle, causing a large error both vertically
%   and horizontally.

%   The markers are at constant opacity, to give an idea of density, and of a marker size that is
%   proportional to their relative intensity. 

opacity = 0.6;
 
figNum = figNum + 1; 
fig1 = figure(figNum);  
set(fig1,'Position',[100,0,1820,1080])
subplot(1,3,1);
spotPositionCellLength(Lcfp, Pcfp, Icfp, 'CFP', 'c', opacity, [floor(leftLim) ceil(rightLim)]);
set(gca,'Color',[0. 0. 0.]);

subplot(1,3,2);
spotPositionCellLength(Lyfp, Pyfp, Iyfp, 'YFP', 'y', opacity, [floor(leftLim) ceil(rightLim)]);
set(gca,'Color',[0. 0. 0.]);

subplot(1,3,3);
spotPositionCellLength(Lrfp, Prfp, Irfp, 'RFP', 'r', opacity, [floor(leftLim) ceil(rightLim)]);
set(gca,'Color',[0. 0. 0.]);

%   This figure shows the spot position versus intensity. THe goal of this is to see if teh
%cell length is correlated to the spot intensity, but the results so far clearly show that this 
%is mostly random. A polynomial fit doesn't tell us much; you need a model before doing so.


figNum = figNum + 1; 
fig2 = figure(figNum);
set(fig2,'Position',[100,0,1820,1080])
subplot(1,3,1) 
spotPositionIntensity(Pcfp, Icfp, 'b', 'CFP');
subplot(1,3,2) 
spotPositionIntensity(Pyfp, Iyfp, 'y', 'YFP');
subplot(1,3,3) 
spotPositionIntensity(Prfp, Irfp, 'r', 'RFP');
 
 
% This shows the spots again, but with histograms showing the distribution of spots.
%       bin size is calculated by Sturges' formula
figNum = figNum + 1; 
fig3 = figure(figNum);
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


% heat maps and 3d plots

% figure 4 contains heat maps
% figure 5 contains 3d plots

figNum = figNum + 1; 
fig4 = figure(figNum);
set(fig4,'Position',[100,0,1820,1080])

subplot(1,3,1);
heatMap( Lcfp, Pcfp, Icfp, calibration(1), 'CFP');

subplot(1,3,2);
heatMap( Lyfp, Pyfp, Iyfp, calibration(2), 'YFP');

subplot(1,3,3);
heatMap( Lrfp, Prfp, Irfp, calibration(3), 'RFP');


%figure  was just the hist3 of the heatmap. 

%% Full cell intensity vs. celllength


figNum = figNum + 1; 
fig6 = figure(figNum);
set(fig6,'Position',[100,0,1820,1080])

subplot(1,3,1);
plotDataCFP = stoichiometryLength(Lcfp, Fcfp, Icfp, calibration(1), 'CFP', 'TetR');
subplot(1,3,2);
plotDataYFP = stoichiometryLength(Lyfp, Fyfp, Iyfp, calibration(1), 'YFP', 'Tus');
subplot(1,3,3);
plotDataRFP = stoichiometryLength(Lrfp, Frfp, Irfp, calibration(1), 'RFP', 'DnaN');
 
if 0 % no idea what this is supposed to be for
    plotyfp(1,:) = Lyfp;
    plotyfp(2,:) = Fyfp;
    plotyfp(3,:) = Iyfp;

    plotyfp = unique(plotyfp','rows')';

    [N_full,Edges_full,mid_full,loc_full]=histcn([plotyfp(1,:)' plotyfp(2,:)'],thisedge3{1},thisedge3{2});
    [N_spot,Edges_spot,mid_spot,loc_spot]=histcn([plotyfp(1,:)' plotyfp(3,:)'],thisedge3{1},thisedge3{2});
end
    
%abort( 'End of figures that seemed correct. -%s', 'Josko');