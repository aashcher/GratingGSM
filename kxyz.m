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
%% wavevector projections and moduli in a given pair of media:
function [kz1, kz2, kx, ky, kxy] = kxyz(no, k0, kg, eps)
	[kx,ky] = meshgrid(k0(1) + kg(1)*(linspace(1,no(1),no(1)) - ceil(no(1)/2)), ...
										 k0(2) + kg(2)*(linspace(1,no(2),no(2)) - ceil(no(2)/2)));
	kx = (reshape(kx,1,[])).';
	ky = (reshape(ky,1,[])).';
	kxy = kx.^2 + ky.^2;

	kz1 = sqrt(eps(1) - kxy);
	kz2 = sqrt(eps(2) - kxy);
	ind = angle(kz1) < -1e-12;
	kz1(ind) = -kz1(ind);
	ind = angle(kz2) < -1e-12;
	kz2(ind) = -kz2(ind);

	kxy = sqrt(kxy);
end