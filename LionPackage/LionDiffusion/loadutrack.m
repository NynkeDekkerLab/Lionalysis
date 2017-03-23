function [A] = loadutrack(N, dataDirectory)
    %LOADUTRACK Summary of this function goes here
    %   Detailed explanation goes here

    if nargin<1
        N=[];
    end

    %mainstring = sprintf('%s\utrack\utrackResults\', dataDirectory);
    %followup='\TrackingPackage\tracks\Channel_1_tracking_result.mat'; 

    for i=1:length(N)
        
        %Nstr=num2str(N(i));
        %tracksFinal{i}=load(strcat(mainstring,Nstr,followup));
        fprintf('File number %d \n ', N(i));
        fileName = sprintf('%s%d/utrackResults/TrackingPackage/tracks/Channel_1_tracking_result.mat', dataDirectory, N(i))
        tracksFinal{i} = load( fileName );
    end

    A=tracksFinal;
end

