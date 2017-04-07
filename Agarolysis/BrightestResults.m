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

%reformat data
for i=1:Nexp 
    Ncells{i}=size(E{i}.DataStruct,2); 
    for j=1:Ncells{i}    
        fprintf('Experiment %d, Cell %d, cellTotalNumber %d: \n', i,j, cellTotalNumber);
        for k=1:3
            channel = 'CFP';
            if k == 2
                channel = 'YFP';
            end
            if k == 3
                channel = 'RFP' ;
            end  

            fprintf('\tChannel %s:\n', channel);

            spots = E{i}.DataStruct(k,j).ld;
            fprintf('\t\tFound %d spots.\n', length(spots));


            subplot(1,3, k); 
            colormap hot;
            view([0 90])
        end  
        pause; 
        cellTotalNumber = cellTotalNumber + 1; 
    end
end
close all;
fprintf('Data has been loaded and re-formatted. Next, figures.\n'); 