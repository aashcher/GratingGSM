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
%% get index in the amplitude or efficiecy vector from diffraction orders and the polarization
function ind = index(ix, iy, no, pol)
	if strcmp(pol, 'TE')
		ind = (ceil(no(1)/2)-1+ix)*no(2) + ceil(no(2)/2) + iy;
	elseif strcmp(pol, 'TM')
		ind = no(1)*no(2) + (ceil(no(1)/2)-1+ix)*no(2) + ceil(no(2)/2) + iy;
	else
		error('index error');
end