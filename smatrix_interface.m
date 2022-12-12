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
%% interface S-matrix:
function [SM] = smatrix_interface(kz1, kz2, eps1, eps2, polarization)
	no = numel(kz1);

	if strcmp(polarization, 'TE')
		SM = zeros(no,2,2);

		SM(:,1,1) = (kz1-kz2)./(kz1+kz2);
		SM(:,2,1) = 1 + SM(:,1,1);
		SM(:,2,2) = -SM(:,1,1);
		SM(:,1,2) = 1 + SM(:,2,2);
	elseif strcmp(polarization, 'TM')
		SM = zeros(no,2,2);

		SM(:,1,1) = (eps2*kz1-eps1*kz2)./(eps2*kz1+eps1*kz2);
		SM(:,2,1) = 1 + SM(:,1,1);
		SM(:,2,2) = -SM(:,1,1);
		SM(:,1,2) = 1 + SM(:,2,2);
	elseif strcmp(polarization, 'TEM')
		ind_e = 1:no;
		ind_h = no+1:2*no;
		SM = zeros(2*no,2,2);
			% TE:
		SM(ind_e,1,1) = (kz1-kz2)./(kz1+kz2);
		SM(ind_e,2,1) = 1 + SM(ind_e,1,1);
		SM(ind_e,2,2) = -SM(ind_e,1,1);
		SM(ind_e,1,2) = 1 + SM(ind_e,2,2);
			% TM:
		SM(ind_h,1,1) = (eps2*kz1-eps1*kz2)./(eps2*kz1+eps1*kz2);
		SM(ind_h,2,1) = 1 + SM(ind_h,1,1);
		SM(ind_h,2,2) = -SM(ind_h,1,1);
		SM(ind_h,1,2) = 1 + SM(ind_h,2,2);
	end
end


