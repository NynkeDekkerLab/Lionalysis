function [] = spotPositionCount(p, I, colour, channel)


	numBins = log(length(p))/log(2)+1;

	edgeLocations = (0:numBins)/numBins;

	[numbin,edges] = histcounts(p,edgeLocations);
	normal = max(numbin)/80;

	binCenter = diff(edges);
	binCenter = cumsum(binCenter) - binCenter(1)/2;

	hold on
		bar(binCenter, numbin/normal, 'w')
		%scatter(p,I,'b','x');
		opacity = 0.5;
		scatter(p,I,20.0, colour,'filled', ...
			'MarkerFaceAlpha', opacity, ...
			'MarkerEdgeAlpha', opacity);
		%plot(binCenter,numbin/normal,'k','LineWidth',3)
	hold off

	xlabel('Position in cell');
	ylabel('Amount of spots (normalized)');
	title(sprintf('Agar data: %s', channel)) ;
	axis([0 1 -0.1 ceil(max(I)/10)*10]);
	set(gca,'FontSize',16);

end