function [] = heatMap( l, p, I, calibration, channel)

	%Filter on Intensity
	I=I*calibration;
	[bin_cfp,idx_cfp]=find(I>calibration*2); 

	pixelUnit = 0.159;

	num(1,:) = l*pixelUnit; 
	num(2,:) = p;

	num(:,idx_cfp)=num(:,idx_cfp);
	numZ(1,:)=nonzeros(num(1,:));
	numZ(2,:)=nonzeros(num(2,:));

	filtedSIgnal=size(numZ,2)/size(I,2); % Percentage of total signal
  
	numHorizontalBins = sturgesFormula(l);
	numVerticalBins = sturgesFormula(p);

	binArray{1} = linspace(15,35,numHorizontalBins+1)*pixelUnit;


	binArray{2} = (0:numVerticalBins)/numVerticalBins;

	intensityMap = hist3(numZ','Edges',binArray);
	intensityMap = intensityMap';
	maximumIntensity = max(intensityMap);

	for i=1:size(intensityMap,2)
		intensityMap(:,i)=intensityMap(:,i)./maximumIntensity(i);
	end
	pcolor(binArray{1},(binArray{2}),intensityMap);
	xlabel('Cell Length');
	ylabel('Position in Cell');
	title(sprintf('Agar data: %s', channel));
	grid on
	set(gca,'FontSize',16)  
    colormap hot;

end