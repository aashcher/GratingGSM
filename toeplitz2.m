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
%% return a 2D Toeplitz matrix
function [T] = toeplitz2(M, nx, ny)
	CT = cell(1,2*nx-1);
	ind = toeplitz(linspace(nx,2*nx-1,nx),flip(linspace(1,nx,nx)));		
	for i = 1:2*nx-1
		CT{1,i} = toeplitz(M(ny:2*ny-1,i),M(ny:-1:1,i));
	end
	T = cell2mat(CT(ind));
end