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
%% calculate 2D Fourier coefficients for the permittivity distribution of multicylinder grating:
function FE = feps_cyl(no, rpx, rpy, cpx, cpy, eps_c, eps_m, eps_b)
	nc = numel(rpx);
	if nc ~= numel(rpy) || nc ~= numel(cpx) || nc ~= numel(cpy) || nc ~= numel(eps_c)
		error('incorrect parameters in function calc_emntd_ncyl');
	end

	FE = cellmat(2,1,2*no(2)-1,2*no(1)-1);
	
	ix = linspace(1,no(1)-1,no(1)-1);
	iy = linspace(1,no(2)-1,no(2)-1);
	[IX,IY] = meshgrid(ix,iy);

	for i = 1 : nc
		fx = rpy(i)*besselj(1,(2*pi*rpx(i))*ix) ./ ix;
		fy = rpx(i)*besselj(1,(2*pi*rpy(i))*iy) ./ iy;
		expx = exp(-(1i*2*pi*cpx(i))*ix);
		expy = exp(-(1i*2*pi*cpy(i))*iy);
		FXY = sqrt((rpx(i)*IX).^2 + (rpy(i)*IY).^2);
		FXY = (rpx(i)*rpy(i)) * besselj(1, (2*pi)*FXY) ./ FXY;
		EXPX = exp(-(1i*2*pi*cpx(i))*IX);
		EXPY = exp(-(1i*2*pi*cpy(i))*IY);

		M = zeros(2*no(2)-1,2*no(1)-1);

		M(no(2),no(1)) = pi*rpx(i)*rpy(i);

		M(no(2)+1:2*no(2)-1,no(1)) = fy .* expy;
		M(no(2)-1:-1:1,no(1)) = fy .* conj(expy);
		M(no(2),no(1)+1:2*no(1)-1) = fx .* expx;
		M(no(2),no(1)-1:-1:1) = fx .* conj(expx);

		M(no(2)+1:2*no(2)-1,no(1)+1:2*no(1)-1) = FXY .* EXPX .* EXPY;
		M(no(2)+1:2*no(2)-1,no(1)-1:-1:1) = FXY .* conj(EXPX) .* EXPY;
		M(no(2)-1:-1:1,no(1)+1:2*no(1)-1) = FXY .* EXPX .* conj(EXPY);
		M(no(2)-1:-1:1,no(1)-1:-1:1) = FXY .* conj(EXPX .* EXPY);

		FE{1} = FE{1} + ((eps_c(i) - eps_m)/eps_b) * M;
		FE{2} = FE{2} + (eps_b/eps_c(i) - eps_b/eps_m) * M;
	end
	
	FE{1}(no(2),no(1)) = FE{1}(no(2),no(1)) + eps_m/eps_b;
	FE{2}(no(2),no(1)) = FE{2}(no(2),no(1)) + eps_b/eps_m;
end
