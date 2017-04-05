function [] = spotPositionIntensity(p,  I,  colour,  channel) 
	myfit = polyfit(p, I, 4);
	x = 0:0.001:1;
	y = polyval(myfit, x);
	% It might be more significant to look at P(I | p), the distribution of intensity
	%	given a certain cell length. But the results seem to clearly indicate randomness.
	hold on
	    scatter(p, I, colour, 'x');
	    plot(x, y, 'k', 'LineWidth', 3)
	hold off

	xlabel('Position in cell');
	ylabel('Spot Intensity'); 
	title(sprintf('Agar data: %s', channel));
	axis([0 1 -0.1 ceil(max(I)/10)*10])
	set(gca, 'FontSize', 16)

end