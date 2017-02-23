%% Define experiment and user
% Each seperate function wants to have these two numbers, which it uses to
% retrieve settings. This seems quite convoluted but appears to have the
% intention of allowing users to stop during analysis and move on at a
% later point. Note that the way this works is by calling
% Axxx_Set_experiment. This accepts the user and exp, and the user accepts
% the exp. 

% The scripts mkae use of MIJ, which is a matlab API for FIJI. FIJI is what ImageJ wraps about.
% As a result, numerous java scripts are called through MIJi.m
exp ='test';
user = 'Josko';

if user == 'Josko'
	if JoskoInit ~= 1
	    %Start in git repository
	    cd('C:\GitHub\KymoCode');
	    % Add the paths of the supporting functions.
	    addpath('SupportingFunctions');
	    addpath('ExperimentSpecificSettings'); 
	    % Load and start DIP
	    addpath('C:\Program Files\DIPimage 2.8');
	    dipstart; 
    end
end;
fprintf('Starting A01: Spine..\n');
%% Section in Matlab. The goal of this particular section is printed on the next line.
fprintf('Loading A010: Prepare Images.\n Fluorescence Analysis requires background noise reduction.\n A010 applies a ten-point rolling background substraction.\n Furthermore, it handles image translation, rotation and scaling.\n These are necessary for proper alignment of the images from Brightfield and Fluorescence.\n'); 
% Sets experiment
% Rolling Ball correction. Rolling ball corrections attempt to correct for uneven illumination. It takes te average intensity of a radial surface around the pixel, and substract it.
% It is a standard plugin for ImageJ.
A010_Prepare_Images_After_NIS(user,exp);

fprintf('End of section.\n\n\n');
%% Section in Matlab. The goal of this particular section is printed on the next line.
fprintf('Loading A20: Save images to Workspace.\n This stores the images (both original and corrected) from both the brightfield and fluorescent cameras in a single tif stack.\n It then loads the images into Matlab memory, which is beneficial for computation times.\n');
% Reads Tiff (stacks) images and puts the data on the workspace.
A020_Images_SaveImageWorkSpace(user,exp);

fprintf('End of section.\n\n\n');
%% Section in Matlab. The goal of this particular section is printed on the next line.
fprintf('Loading A30: Drift calculation.\n In time-lapse imaging, the spatial drift of the smaple might be significant due to the thermal gradients and other factors.\n Due tot the limitiations of teh temperature control system in the experimental setup, it is subject to some thermal expansion and contraction.\n A030 finds the drift vector of the sample from its initial position of the course of the measurement.\n Usually, the most significant sourfce of drift is translational motion, which means that all points in the image share the same drift regardless of their position in the field of view.\n The localisation event-based cross-correlation methods exploited in A030 utilize the cross-correlation function to estimate the drift.\n The maximal value of teh cross correlation function corresponds exactly to the extent of the sample drift between two time points.\n The drift vector is estimated purely from bright field imaging.\n');
%
A030_Images_FindDriftVector(user,exp)
close all

fprintf('End of section.\n\n\n');
%% Section in Matlab. The goal of this particular section is printed on the next line.
fprintf('Loading A040: Kymomaker.\n A040 needs to know the exact position of a growth channel, determined through user input.\n The user selects a singular growth channel and the other channels are then estimated using the properties of hte microfluidic device.\n The Kymographs are constructed from the Fluorescent images for DNAN.\n');

initval = A040_Images_KymoMaker(user,exp);
WorkspaceOutName = initval.WorkspaceOutName;
close all

fprintf('End of section.\n\n\n');
%% Section in Matlab. The goal of this particular section is printed on the next line.
fprintf('Loading A060: Replication Cluster.\n This module identifies the loading and unloading events of DnaN in time and space via a series of user inputs.\n A060 prompts the user to click on a set of points on the Kymograph.\n First, the user clicks on the localization point of a DnaN spot, which first becomes visible.\n Thereafter, the user clicks to identify the delocalization point of DnaN in hte same cell, where the fluorescence spot disappears (indicative of replication termination).\n Finally, the user is prompted to indicate the two new starting points of DnaN activity in the daughter cells.\n ');

if exist('WorkspaceOutName','var')
    A060_RepliCluster00_Click(user,exp,WorkspaceOutName)
else
    A060_RepliCluster00_Click(user,exp)
end
close all

fprintf('End of section.\n\n\n');
%% Section in Matlab. The goal of this particular section is printed on the next line.
fprintf('A100 creates Bacpics using the replication cycles defined in A060.\n However, it consists of several submodules.\n First, center of mass tracking is performed on the Kymographs to find the traces of the DnaN spots on them.\n The localisation and delocalisation of each spot were indicated by the user in A060 and are used as boundary conditions.\n The traces are stored together with the original Kymographs.\n Next, the timepoints of cell birth and division associated with their replication cycle are detected from the brightfield image.\n The cell edges are detected for each frame.\n Then, a more detailed analysis is made on the precise times of replication initiation and termination as opposed to the manual clicks, based on setep-fitting the spot signals.\n ');

if exist('WorkspaceOutName','var')
    A100_Processing_ReplicationAutoShell(user,exp,WorkspaceOutName)
else
    A100_Processing_ReplicationAutoShell(user,exp)
end

fprintf('End of section.\n\n\n');
%% Section in Matlab. The goal of this particular section is printed on the next line.
fprintf('Save Bacpics to workspace .\n');

[Bacpics, MDchans, FLchans, Cells] = Get_Bacpics(user,exp);
fprintf('End of section.\n');
%% Section in Matlab. The goal of this particular section is printed on the next li\n\nne.
fprintf('\n Bacpics saved to Workspace');

fprintf('End of section.\n\n\n');
%% Section in Matlab. The goal of this particular section is printed on the next line.
fprintf('.\n');

initval=A001_Images_Set_Experiment(user,exp);
finddif = strcmp(initval.difchan,initval.viewchan);
viewFLchans = {initval.viewchan{finddif},initval.viewchan{~finddif}};
IPTP = 1;

for chan = 1:FLchans;
    IPTPvalue(chan) = LionfitcontrolUI(1,Bacpics{1,chan},viewFLchans,Cells(1),chan,'Kymo');
end


fprintf('End of section.\n\n\n');
%% Section in Matlab. The goal of this particular section is printed on the next line.
fprintf('Lionfit .\n');

initval.Lionpath = strcat(initval.Kymopath,initval.OSslash,'LionFit',initval.OSslash);
addpath(genpath(strcat(initval.Lionpath,'150917V')));
addpath(genpath(strcat(initval.Lionpath,'gaussmlev2')));

if exist('IPTPvalue','var')
    GaussFitSimedit(user,exp,IPTPvalue)
else
    GaussFitSimedit(user,exp)
end
