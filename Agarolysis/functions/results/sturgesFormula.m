function [numBins] = sturgesFormula(p)
	if 0
		numBins = log(length(p))/log(2)+1;
	else
		warning('AgarResults::sturgesFormula$:Sturges formula is insufficient. Using square root law.');
		numBins = ceil( sqrt(length(p)));
	end
end