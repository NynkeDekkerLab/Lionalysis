if 0
    experimentNames = [1 2 3 4 5 6]; %experiment names
    LionStart(experimentNames);  % importeren trajectories which utrack detects, saves in data/Results
    %you caan change PAMCherry as a name in the General_set used for Vbspt   
    Dchoplot('General_set'); %probability diff constants


    PlotTrajs; %asks for chopped trajs or not, set threshold between red/blue coloured
PlotDT; %dwell time analysis
end

cd('DiffDistr');
MainScript