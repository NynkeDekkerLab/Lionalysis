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
celllength = []; 

%   pure loading
j=1; 
for i=exps;
    E{j}=load(strcat(folder,slash,num2str(i),slash,'Results.mat')); 
    j=j+1;
end

allCFP_L = [];

Nexp=size(E,2);


cellTotalNumber = 0;

twoDimGauss = @(x,y,mux, stdx, muy, stdy) 1 / 4 / pi / stdx / stdy .* exp(-(x-mux).^2/4/stdx^2 - (y-muy).^2/4/stdy^2);
figNumber = 1;

inspectFigure = figure(figNumber);
set(inspectFigure,'Position',[0,0,1680,1080]) 

inspectOpacity = 0.7;
%spotPositionCellLength(Lcfp, Pcfp, Icfp, 'CFP', 'c', opacity, [floor(leftLim) ceil(rightLim)]);

cellLengthCFP = [];
cellLengthYFP = [];
cellLengthRFP = [];

positionCFP = [];
positionYFP = [];
positionRFP = [];

intensityCFP = [];
intensityYFP = [];
intensityRFP = [];

button = questdlg('Display each cell?','Settings', 'No');

displayCells = strcmp( button, 'Yes');


%for each experiment
for i=1:Nexp 
    Ncells{i}=size(E{i}.DataStruct,2); 

    for j=1:Ncells{i}    
        fprintf('Experiment %d, Cell %d, cellTotalNumber %d: \n', i,j, cellTotalNumber);
        clf;

        cyanSpots = [];
        yellowSpots = [];
        redSpots = [];

        cyanIntensity = [];
        yellowIntensity = [];
        redIntensity = [];

        %for each colour channel
        for k=1:3
            channel = 'CFP';
            if k == 2
                channel = 'YFP';
            end
            if k == 3
                channel = 'RFP' ;
            end  

            spotIntense = [0.00];

            fprintf('\tChannel %s:\n', channel);

            cellLength = E{i}.DataStruct(k,j).CellLength;

            spots = E{i}.DataStruct(k,j).ld;
            fprintf('\t\tFound %d spots.\n', length(spots)); 

            for spotNum=1:length(spots)
                spot = spots(spotNum);
                spot = spot{:}; 

                if spot(1) > spotIntense(1)
                    spotIntense = spot;
                end


                if k == 1
                    cyanSpots = [cyanSpots, [spot(2) spot(4)]'];
                    cyanIntensity = [cyanIntensity, spot(1)];
                elseif k==2
                    yellowSpots = [yellowSpots, [spot(2) spot(4)]'];
                    yellowIntensity = [yellowIntensity, spot(1)];
                else 
                    redSpots = [redSpots, [spot(2) spot(4)]'];
                    redIntensity = [redIntensity, spot(1)];
                end 
            end  
            if k == 1
                cellLengthCFP = [cellLengthCFP, cellLength];
                intensityCFP = [intensityCFP, spot(1)];
                positionCFP = [positionCFP, spot(2)];
            elseif k==2
                cellLengthYFP = [cellLengthYFP, cellLength];
                intensityYFP = [intensityYFP, spot(1)];
                positionYFP = [positionYFP, spot(2)];
            else 
                cellLengthRFP = [cellLengthRFP, cellLength];
                intensityRFP = [intensityRFP, spot(1)];
                positionRFP = [positionRFP, spot(2)];
            end 
        end  

        minNumSpots = 1e5;

        if length(cyanSpots) < minNumSpots
            minNumSpots = length(cyanSpots);
        end

        if length(yellowSpots) < minNumSpots
            minNumSpots = length(yellowSpots);
        end

        if length(redSpots) < minNumSpots
            minNumSpots = length(redSpots);
        end

        if minNumSpots > 0
            cyanSpots = cyanSpots ./ E{i}.DataStruct(3,j).CellLength;
            yellowSpots = yellowSpots ./ E{i}.DataStruct(3,j).CellLength;
            redSpots = redSpots ./ E{i}.DataStruct(3,j).CellLength;

            cyanIntensity = 400*cyanIntensity / max(cyanIntensity);
            yellowIntensity = 400*yellowIntensity / max(yellowIntensity);
            redIntensity = 400*redIntensity / max(redIntensity);

            if displayCells
                hold on;
                    scatter(cyanSpots(2,:), cyanSpots(1,:), cyanIntensity, 'c', ...
                        'filled', ...
                        'MarkerFaceAlpha', inspectOpacity, ...
                        'MarkerEdgeAlpha', inspectOpacity);
                    scatter(yellowSpots(2,:), yellowSpots(1,:), yellowIntensity, 'y', ...
                        'filled', ...
                        'MarkerFaceAlpha', inspectOpacity, ...
                        'MarkerEdgeAlpha', inspectOpacity);
                    scatter(redSpots(2,:), redSpots(1,:), redIntensity, 'r', ...
                        'filled', ...
                        'MarkerFaceAlpha', inspectOpacity, ...
                        'MarkerEdgeAlpha', inspectOpacity);
                hold off;

                
                xlim([-.2 .2]);
                ylim([0 1.0]);
                set(gca,'Color',[0. 0. 0.]);

                pause(.3); 
            end
        else
            fprintf('Skipping data points because one or more channels detected zeor spots.\n');
        end
        cellTotalNumber = cellTotalNumber + 1; 
    end
end
close all;
fprintf('Data has been loaded and re-formatted. Next, figures.\n'); 
 
figNumber = figNumber + 1; 
fig0 = figure(figNumber);

lengthArray = [cellLengthCFP, cellLengthYFP, cellLengthRFP];

histogram( lengthArray, linspace( floor(min(lengthArray)), ceil(max(lengthArray)), log(length(lengthArray))/log(2)+1));

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
 
figNumber = figNumber + 1; 
fig1 = figure(figNumber);  
set(fig1,'Position',[100,0,1820,1080])

intensityCFP = 400 * intensityCFP / max(intensityCFP);
intensityYFP = 400 * intensityYFP / max(intensityYFP);
intensityRFP = 400 * intensityRFP / max(intensityRFP);

normalisedPositionCFP = positionCFP./cellLengthCFP;
normalisedPositionYFP = positionYFP./cellLengthYFP;
normalisedPositionRFP = positionRFP./cellLengthRFP;

subplot(1,3,1);
spotPositionCellLength(cellLengthCFP, normalisedPositionCFP, intensityCFP, 'CFP', 'c', opacity, [floor(leftLim) ceil(rightLim)]);
set(gca,'Color',[0. 0. 0.]); 

subplot(1,3,2);
spotPositionCellLength(cellLengthYFP, normalisedPositionYFP, intensityYFP, 'YFP', 'y', opacity, [floor(leftLim) ceil(rightLim)]);
set(gca,'Color',[0. 0. 0.]); 

subplot(1,3,3);
spotPositionCellLength(cellLengthRFP, normalisedPositionRFP, intensityRFP, 'RFP', 'r', opacity, [floor(leftLim) ceil(rightLim)]);
set(gca,'Color',[0. 0. 0.]); 


figNumber = figNumber + 1; 
fig1 = figure(figNumber);  
set(fig1,'Position',[100,0,1820,1080])

subplot(1,3,1);
spotDistanceCellLength(cellLengthCFP, normalisedPositionCFP,cellLengthRFP, normalisedPositionRFP, 'CFP - RFP', 'c', opacity, [floor(leftLim) ceil(rightLim)]);
set(gca,'Color',[0. 0. 0.]); 

subplot(1,3,2);
%spotDistanceCellLength(cellLengthYFP, normalisedPositionYFP,cellLengthRFP, normalisedPositionRFP, 'YFP - RFP', 'r', opacity, [floor(leftLim) ceil(rightLim)]);
spotDistanceCellLength(cellLengthCFP, normalisedPositionCFP,cellLengthRFP, rand(1,length(normalisedPositionRFP)), 'CFP - Rand', 'r', opacity, [floor(leftLim) ceil(rightLim)]);
set(gca,'Color',[0. 0. 0.]);

subplot(1,3,3);
spotDistanceCellLength(cellLengthYFP, normalisedPositionYFP,cellLengthCFP, normalisedPositionCFP, 'CFP - YFP', 'g', opacity, [floor(leftLim) ceil(rightLim)]);
set(gca,'Color',[0. 0. 0.]); 