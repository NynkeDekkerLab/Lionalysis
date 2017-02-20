function [] = AgarCD(path)
% CDs to a path without throwing errors and breaking everything.
	try
		cd(path)
	catch
		fprintf('Failed to CD to %s. \n', path);
	end
end