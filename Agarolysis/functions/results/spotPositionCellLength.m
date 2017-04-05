function [] = spotPositionCellLength(Lcfp, Pcfp, Icfp, channel, colour) 
	if 0
		hold on
		scatter(single(Lcfp),Pcfp,Icfp,colour,'filled');
		myfit=polyfit(Lcfp,Pcfp,4);
		x=15:0.1:45;
		y=polyval(myfit,x);
		plot(x,y,'r','LineWidth',5)
		legend('spot scatter', 'polynomial fit P_4')
		hold off
	else
		scatter(single(Lcfp),Pcfp,Icfp,colour,'filled');
	end

	xlabel('Cell Length');
	ylabel('Normalized spot position in cell'); 
	title(sprintf('Agar data: %s', channel))
	axis([12 43 -0.1 1.1])
	set(gca,'FontSize',16)	
end