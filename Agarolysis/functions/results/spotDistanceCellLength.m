function [] = spotDistanceCellLength(firstLength, firstPosition, secondLength, secondPosition, channel, colour, opacity, horizontalLimit) 
% scatter plot for cell length versus P1 P2 distances

    meanLength = 0.5 * ( firstLength + secondLength);

    longDistance = abs(secondPosition - firstPosition);


    scatter( meanLength, longDistance, colour,'filled', ...
        'MarkerFaceAlpha', opacity, ...
        'MarkerEdgeAlpha', opacity); 

    xlabel('Cell Length');
    ylabel('Spot distance normalised to cell length'); 
    title(sprintf('Agar data: %s', channel)) 

    left = horizontalLimit(1);
    right = horizontalLimit(2);

    axis([left right -0.1 1.1])
    set(gca,'FontSize',16)  
end