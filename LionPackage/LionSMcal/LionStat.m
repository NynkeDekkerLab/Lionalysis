clc
clear all
close all
%This file displays results.
%% Load Results

TraceNumbers=[9,11];

init.OSslash = '/';

fprintf('Select Final Results Folder');
init.resultspath = uigetdir(pwd,'Select Results Folder');
    
if init.resultspath == 0;
    init.resultspath = '/Users/rleeuw/Work/Data/170111_Tus-SMcal/gain100/FinalResults';
end

init.resultspath = strcat(init.resultspath,init.OSslash);

Results=cell(length(TraceNumbers),1);

w=1;
for i=TraceNumbers
    Results{w}=load(strcat(init.resultspath,'SMCResult',num2str(i)));
    w=w+1;
end

%% Collect all integrated intensity data
IItot=[];

for i=1:size(Results,1)
    IItot=[IItot;Results{i}.SMCResult.TraceVal];
end

%% plot histogram
Nbins=.5*( sqrt(length(IItot)) + log(length(IItot))/log(2) + 1);
Xrightbound=3000;
IItotc=histc(IItot,Nbins);
fig1=figure(1);
normalisedI=IItotc/(sum(IItotc));
bar(bins, normalisedI)

meanValue = 0;
for i = 1:length(normalisedI)
    meanValue = meanValue + (bins(i)+bins(i))/2 * normalisedI(i);
end
title(sprintf('Mean value  %.3f', meanValue));
set(gca,'fontsize',18)
axis([0 Xrightbound 0 0.2])
xlabel('Integrated Intensity (counts)')
ylabel('Probability Density (1/counts)')
