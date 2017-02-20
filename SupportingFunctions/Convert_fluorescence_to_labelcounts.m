function correctedlabels=Convert_fluorescence_to_labelcounts(countsperlabel,fluo,frames);
%function corrects for bleaching; JacobKers 2013

countsperbetaclamp=countsperlabel;
%This is the conversion unit when no bleaching is present

decreasefactor=2.1;  
%This is a factor describing how much the fluorescence
%of bacteria will have dropped when bleaching and expression balannce in a 'steady-state' signal 

decay=30; 
%This is the time constant of bleaching in frame units, assuming
%exponential decay to the steady-state level

%Above factors must be determined elsewhere.

if nargin<2  %if TEST MODE
    close all;
    frames=[10:60];
    labels=0*frames+100;   %constant 100 labels
    bleachfactorperframe=((decreasefactor-1)*exp(-frames/decay)+1)/decreasefactor;
    bleachfactorperframe=((decreasefactor-1)*exp(-frames/decay)+1)/decreasefactor;
    fluo=countsperbetaclamp*labels.* bleachfactorperframe;     
end
%-------------------------------------------------------------------

correctfactor=1./(((decreasefactor-1)*exp(-frames/decay)+1)/decreasefactor)';
%correctfactor=2;
uncorrectedlabels=fluo/countsperbetaclamp;

%correctfactor
%uncorrectedlabels
correctedlabels=correctfactor.*uncorrectedlabels;

if nargin<2 %TEST MODE
     plot(frames,labels,'-'); hold on;
     plot(frames,uncorrectedlabels,'-o');
     plot(frames,correctedlabels,'ro');
     legend('original labels', 'calculated labels, no bleachcorrection', 'calculated labels, bleachcorrected' );
end