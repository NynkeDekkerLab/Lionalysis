function [] = spotPositionCellLength(L, P, I, channel, colour, opacity, horizontalLimit) 
	if 0
		hold on
		scatter(single(L),P,I,colour,'filled');
		myfit=polyfit(L,P,4);
		x=15:0.1:45;
		y=polyval(myfit,x);
		plot(x,y,'r','LineWidth',5)
		legend('spot scatter', 'polynomial fit P_4')
		hold off
	else
		scatter(single(L),P,I/mean(I)*40.0, colour,'filled', ...
			'MarkerFaceAlpha', opacity, ...
			'MarkerEdgeAlpha', opacity);
	end

	xlabel('Cell Length');
	ylabel('Normalized spot position in cell'); 
	title(sprintf('Agar data: %s', channel))

	%left = mean(L) - std(L);
	%right = mean(L) + std(L);

	left = horizontalLimit(1);
	right = horizontalLimit(2);

	axis([left right -0.1 1.1])
	set(gca,'FontSize',16)	
end