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
%%
function [C] = circulant2(M, nx, ny)
	CC = cell(1,2*nx-1);
	for i = 1:nx
		CC{i} = [M(ny:-1:1, nx+1-i); M(2*ny-1:-1:ny+1, nx+1-i)];
	end
	for i = 1:nx-1
		CC{nx+1+i} = [M(ny:-1:1, 2*nx-i); M(2*ny-1:-1:ny+1, 2*nx-i)];
	end
	C = cell2mat(CC);
end
