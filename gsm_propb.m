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
%% calculate plane wave propagation from slices to grating boundaries
function U = gsm_propb(V, method)
	no = method.no(1)*method.no(2);
	ns = method.nsl;
	kh = method.kh;
	zs = method.zs;
	kz = method.kz;
	kx = method.kx;
	ky = method.ky;
	kxy = method.kxy;
	zc = method.zs;
	epsb = method.eps_b;
	SM1 = method.SMsubb;
	SM2 = method.SMbsup;
	
	ind1 = 1:no;
	ind2 = no+1:2*no;
	ind3 = 2*no+1:3*no;
	ind4 = 3*no+1:4*no;
	indp = [ind1,ind2];
	indm = [ind3,ind4];

	V = reshape(V, [3*no,ns]);
	
	A = zeros(4*no, ns);
	Ab = zeros(4*no, 1);
	U = zeros(2*no, 2);
	texp = zeros(2*no, ns);
		% x,y,z components to wave amplitudes
	A(ind2,:) = -V(ind3,:) .* kxy ./ kz;
	A(ind4,:) = A(ind2,:);
	A(ind1,:) = (kx .* V(ind1,:) + ky .* V(ind2,:)) ./ kxy;
	A(ind2,:) = A(ind2,:) + A(ind1,:);
	A(ind4,:) = A(ind4,:) - A(ind1,:);
	A(ind1,:) = (-ky .* V(ind1,:) + kx .* V(ind2,:)) ./ (kxy .* kz);
	A(ind3,:) = A(ind1,:);
		% sum amplitudes at layer boundaries
	texp(ind1,:) = exp(1i * kz .* (0.5*kh - zc));
	texp(ind2,:) = texp(ind1,:);
	Ab(indp,1) = sum(A(indp,:) .* texp, 2);
	texp(ind1,:) = exp(1i * kz .* (0.5*kh + zc));
	texp(ind2,:) = texp(ind1,:);
	Ab(indm,1) = sum(A(indm,:) .* texp, 2);
		% self-consistent amplitudes at outer boundaries
	texp(ind1,1) = exp((1i*kh) * kz);
	texp(ind2,1) = texp(ind1,1);
	D = 1 ./ (1 - texp(:,1) .* texp(:,1) .* SM1(:,2,2).*SM2(:,1,1));
	U(:,1) = D .* SM1(:,1,2) .* (Ab(indp,1) .* SM2(:,1,1) .* texp(:,1) + Ab(indm,1));
	U(:,2) = D .* SM2(:,2,1) .* (Ab(indp,1) + Ab(indm,1) .* SM1(:,2,2) .* texp(:,1));
end

%%% END OF FILE %%%