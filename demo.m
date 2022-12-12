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
clc;
format long;
%%
	% fill the method structure
gsm_method.no = [20, 20]; % number of Fourier harmonics
gsm_method.nsl = 30; % number of slices
gsm_method.eps_b = 1; % basis permittivity of the grating layer
gsm_method.restart = []; % GMRes parameter
gsm_method.maxit = 500; % maximum number of iterations for the GMRes
gsm_method.tol = 1e-6; % tolerance for the GMRes

NO = gsm_method.no(1) * gsm_method.no(2); % total number of harmonics
ic = (ceil(gsm_method.no(1)/2)-1) * ...
		gsm_method.no(2)+ceil(gsm_method.no(2)/2); % central harmonic

	% grating containing multiple cylindrical pitches on a period
ncyl = [5, 5]; % numbers of cylinders along X and Y dimensions
Ncyl = ncyl(1) * ncyl(2); % total number of cylinders
	% coordinates of centers relative to the period:
cx = -0.5 + 0.5/ncyl(1) + (1/ncyl(1))*[0:(ncyl(1)-1)];
cy = -0.5 + 0.5/ncyl(2) + (1/ncyl(2))*[0:(ncyl(2)-1)];
[CX, CY] = meshgrid(cx, cy);

grating.layer_type = "2D grating";
grating.type = "cylinder";
grating.period = 0.65 * ncyl;
grating.depth = 0.4;
grating.radius = 0.15 + 0.02*(rand(1, Ncyl)-0.5);
grating.center_x = CX;
grating.center_y = CY;
grating.eps = 3.1^2 * ones(1, Ncyl); % permittivities of cylinders
grating.eps_m = 1; % permittivity of the medium surrounding the cylinders
grating.eps_sub = 1.494^2; % substrate permittivity
grating.eps_sup = 1; % superstrate permittivity

	% incidence angles:
phi_inc = 0;
theta_inc = 10;
	% incidence wavevector projections:
k_inc = [sin(pi*theta_inc/180)*cos(pi*phi_inc/180); sin(pi*theta_inc/180)*sin(pi*phi_inc/180)];
Vinc = zeros(2*NO, 2); % incident amplitude vector
Vinc(ic, 1) = 1; % plane wave incidence from the substrate side

wavelength = 1.55;
wavenumber = 2*pi/wavelength;
kg = [wavelength/grating.period(1), wavelength/grating.period(2)];

incidence.wavelength = wavelength;
incidence.k_inc = k_inc;
incidence.Vinc = Vinc;

	% calculate wavevector projections and moduli for the reciprocal lattice:
[kz, ~, kx, ky, kxy] = kxyz(gsm_method.no, k_inc, kg, [gsm_method.eps_b, gsm_method.eps_b]);
gsm_method.kx = kx;
gsm_method.ky = ky;
gsm_method.kxy = kxy;
gsm_method.kz = kz;
gsm_method.kh = wavenumber * grating.depth;
	% centers of slices along vertical direction:
gsm_method.zs = -0.5*gsm_method.kh + ((1 : gsm_method.nsl)-0.5)*(gsm_method.kh / gsm_method.nsl);
	% integration wights for slices:
gsm_method.zweight = gsm_method.kh / gsm_method.nsl;

	% diagonal S-matrices of substrate and superstrate structure:
[kz1, kz2] = kxyz(gsm_method.no, k_inc, kg, [grating.eps_sub, gsm_method.eps_b]);
SMsubb = smatrix_interface(kz1, kz2, grating.eps_sub, gsm_method.eps_b, 'TEM');
[kz1, kz2] = kxyz(gsm_method.no, k_inc, kg, [gsm_method.eps_b, grating.eps_sup]);
SMbsup = smatrix_interface(kz1, kz2, gsm_method.eps_b, grating.eps_sup, 'TEM');

gsm_method.SMsubb = SMsubb;
gsm_method.SMbsup = SMbsup;

	% 2D Toeplitz vector for grating permittivity Fourier coefficients:
FM = feps_cyl(gsm_method.no, grating.radius/grating.period(1), grating.radius/grating.period(2), ...
							grating.center_x, grating.center_y, grating.eps, grating.eps_m, gsm_method.eps_b);
	% convert to circulant matrices and pre-calculate their FFT:
CM_eps = cell(2,1);
CM_eps{1} = fft2(circulant2(FM{1}, gsm_method.no(1), gsm_method.no(2)));
CM_eps{2} = fft2(circulant2(FM{2}, gsm_method.no(1), gsm_method.no(2)));
gsm_method.CM_eps = CM_eps; % circulant matrix first row

	%% test run:
tic;
[Vout, Veff, balance, rv] = gsm_diffraction(incidence, gsm_method, grating, 1, 0);
toc;
fprintf("balance = %f\n", balance); % accuracy of power conservation
	% display efficiencies:
fprintf("a{0,0} = %f;  a{0,1} = %f\n", Veff(index(0,0,gsm_method.no,'TE'),2), Veff(index(0,1,gsm_method.no,'TE'),2));
fprintf("a{0,-1} = %f;  a{1,0} = %f\n", Veff(index(0,-1,gsm_method.no,'TE'),2), Veff(index(0,1,gsm_method.no,'TE'),2));

return;

% END %