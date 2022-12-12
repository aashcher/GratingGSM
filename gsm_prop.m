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
%% calculate plane wave propagation between slices
function U = gsm_prop(V, method)
	no = method.no(1)*method.no(2);
	ns = method.nsl;
	kh = method.kh;
	zs = method.zs;
	kz = method.kz;
	kx = method.kx;
	ky = method.ky;
	kxy = method.kxy;
	epsb = method.eps_b;
	zc = method.zs;
	SM1 = method.SMsubb;
	SM2 = method.SMbsup;
	
	ind1 = 1:no;
	ind2 = no+1:2*no;
	ind3 = 2*no+1:3*no;
	ind4 = 3*no+1:4*no;
	indp = [ind1,ind2];
	indm = [ind3,ind4];
	
	A = zeros(4*no, ns);
	Ab = zeros(4*no, 1);
	U = zeros(3*no, ns);
		% x,y,z components to wave amplitudes
	A(ind2,:) = -V(ind3,:) .* kxy ./ kz;
	A(ind4,:) = A(ind2,:);
	A(ind1,:) = (kx .* V(ind1,:) + ky .* V(ind2,:)) ./ kxy;
	A(ind2,:) = A(ind2,:) + A(ind1,:);
	A(ind4,:) = A(ind4,:) - A(ind1,:);
	A(ind1,:) = (-ky .* V(ind1,:) + kx .* V(ind2,:)) ./ (kxy .* kz);
	A(ind3,:) = A(ind1,:);
		% sum amplitudes at layer boundaries
	texp = zeros(2*no, ns);
	texp(ind1,:) = exp(1i * kz .* (0.5*kh - zc));
	texp(ind2,:) = texp(ind1,:);
	Ab(indp,1) = sum(A(indp,:) .* texp, 2);
	texp(ind1,:) = exp(1i * kz .* (0.5*kh + zc));
	texp(ind2,:) = texp(ind1,:);
	Ab(indm,1) = sum(A(indm,:) .* texp, 2);
		% self-consistent amplitudes at boundaries
	texp = exp((1i*kh) * kz);
	ta = Ab(ind1,1);
	Ab(ind1,1) = Ab(ind3,1) + Ab(ind1,1) .* SM2(ind1,1,1) .* texp;
	Ab(ind3,1) = ta + Ab(ind3,1) .* SM1(ind1,2,2) .* texp;
	ta = Ab(ind2,1);
	Ab(ind2,1) = Ab(ind4,1) + Ab(ind2,1).*SM2(ind2,1,1) .* texp;
	Ab(ind4,1) = ta + Ab(ind4,1).*SM1(ind2,2,2) .* texp;
	D = 1 ./ (1 - texp .* texp .* SM1(ind1,2,2).*SM2(ind1,1,1));
	Ab(ind1,1) = Ab(ind1,1) .* D .* SM1(ind1,2,2);
	Ab(ind3,1) = Ab(ind3,1) .* D .* SM2(ind1,1,1);
	D = 1 ./ (1 - texp .* texp .* SM1(ind2,2,2).*SM2(ind2,1,1));
	Ab(ind2,1) = Ab(ind2,1) .* D .* SM1(ind2,2,2);
	Ab(ind4,1) = Ab(ind4,1) .* D .* SM2(ind2,1,1);
	texp = exp((1i*(0.5*kh + zc(1))) * kz);
	Ab(ind1,1) = Ab(ind1,1) .* texp;
	Ab(ind2,1) = Ab(ind2,1) .* texp;
	texp = exp((1i*(0.5*kh - zc(end))) * kz);
	Ab(ind3,1) = Ab(ind3,1) .* texp;
	Ab(ind4,1) = Ab(ind4,1) .* texp;
		% amplitudes in layers
	texp = exp((1i*(kh/ns))*kz);
	for i = 1:ns
		ta = A(ind1,i);
		A(ind1,i) = 0.5*A(ind1,i) + Ab(ind1,1);
		Ab(ind1,1) = (Ab(ind1,1) + ta) .* texp;
		ta = A(ind2,i);
		A(ind2,i) = 0.5*A(ind2,i) + Ab(ind2,1);
		Ab(ind2,1) = (Ab(ind2,1) + ta) .* texp;
		ta = A(ind3,ns+1-i);
		A(ind3,ns+1-i) = 0.5*A(ind3,ns+1-i) + Ab(ind3,1);
		Ab(ind3,1) = (Ab(ind3,1) + ta) .* texp;
		ta = A(ind4,ns+1-i);
		A(ind4,ns+1-i) = 0.5*A(ind4,ns+1-i) + Ab(ind4,1);
		Ab(ind4,1) = (Ab(ind4,1) + ta) .* texp;
	end
		% amplitudes to x,y,z components
	ta = (A(ind1,:) + A(ind3,:)) ./ kxy;
	U(ind1,:) = ta .* ky;
	U(ind2,:) = -ta .* kx;
	ta = (A(ind4,:) - A(ind2,:)) .* ((1/epsb) * kz ./ kxy);
	U(ind1,:) = U(ind1,:) + ta .* kx;
	U(ind2,:) = U(ind2,:) + ta .* ky;
	U(ind3,:) = (A(ind2,:) + A(ind4,:)) .* ((1/epsb) * kxy);
	
	U = reshape(U, [3*no*ns,1]);
end

%%% END OF FILE %%%