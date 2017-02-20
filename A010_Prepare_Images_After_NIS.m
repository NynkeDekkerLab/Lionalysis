function A010_Prepare_Images_After_NIS(user,exp)
 
% Sets a number of directories
% Calls an experimental function?
initval = A001_Images_Set_Experiment(user,exp); %define your paths and files

%% Save and process FL images
ColorNum=size(initval.viewchan,2);
for Cidx = 1:ColorNum;
    beamnumber = strcmp(initval.viewchan{Cidx},initval.beampaths(:,1));
    Beampath = initval.beampaths{beamnumber,2};
    RoicotrascaKymo(initval,initval.viewchan{Cidx},Beampath, Cidx)
end

%% Save BF image

BFfolder = strcat(initval.basepath,initval.BFdatapath);
RawBFdir = strcat(BFfolder,initval.rawfolder);
BFtifdir= strcat(BFfolder,initval.BFfiletemplate);

% Add ImageJ and MIJI to path
jarpath = initval.Supportpath;
javaaddpath(strcat(jarpath,'mij.jar'))
javaaddpath(strcat(jarpath,'ij.jar'))

MIJ.start
MIJ.run('Image Sequence...',strcat('open=[',RawBFdir,']'))
MIJ.run('Tiff...',strcat('path=[',BFtifdir,']'))
MIJ.closeAllWindows
MIJ.exit

disp('A10 done')