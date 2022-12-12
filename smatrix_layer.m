%{
Copyright © 2022 Alexey A. Shcherbakov. All rights reserved.
This file is part of the GratingGSM project.

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This sortware is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
%}
%% layer S-matrix:
function [SM] = smatrix_layer(kz, kh, polarization)
	no = numel(kz);

	if strcmp(polarization, 'TE') || strcmp(polarization, 'TM')
		SM = zeros(no,2,2);
		SM(:,2,1) = exp((1i*kh)*(kz.'));
		SM(:,1,2) = SM(:,2,1);
	elseif strcmp(polarization, 'TEM')
		SM = zeros(2*no,2,2);
		SM(1:no,2,1) = exp((1i*kh)*(kz.'));
		SM(1:no,1,2) = SM(1:no,2,1);
		SM(no+1:2*no,2,1) = SM(1:no,2,1);
		SM(no+1:2*no,1,2) = SM(1:no,1,2);
	end
end


