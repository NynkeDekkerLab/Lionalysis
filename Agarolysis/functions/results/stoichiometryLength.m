function [plotData] = stoichiometryLength(L, F, I, calibration, channel, species)
 	plotData(1,:) = L; %length
	plotData(2,:) = F; %full cell intensity
	plotData(3,:) = I*calibration(1); %spot intensity

	plotData = unique(plotData','rows')';

	myfit = polyfit(plotData(1,:),plotData(2,:),1);
	myfit2 = polyfit(plotData(1,:),plotData(3,:),1);

	x=12:0.1:43;
	y=polyval(myfit,x);
	y2=polyval(myfit2,x);

	hold on
		scatter(plotData(1,:),plotData(2,:),'r','o','filled');
		scatter(plotData(1,:),plotData(3,:),'b','o','filled');
		plot(x,y,'g--','LineWidth',1)
		plot(x,y2,'g--','LineWidth',1)
	hold off

	xlabel('Cell Length (px)');
	ylabel('Intensity (-)'); 
	title(sprintf('%s: %s stoichiometry vs. Length', channel, species))
	axis([12 43 -0.1 2e5])
	set(gca,'FontSize',16)
	legend('Full cell intensity','Spot intensity') 
end