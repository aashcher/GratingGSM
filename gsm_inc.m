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
%% calculate amplitudes of the incident field in all slices
function U = gsm_inc(V, method)
	no = method.no(1)*method.no(2); % orders
	ns = method.nsl; % slices
	kh = method.kh; % grating depth
	zs = method.zs; % slice coords, row
	kz = method.kz; % kz(no,1), column
	kx = method.kx;
	ky = method.ky;
	kxy = method.kxy;
	epsb = method.eps_b;
	SM1 = method.SMsubb;
	SM2 = method.SMbsup;

	Va = zeros(4*no, ns);

	te1 = exp((1i*kh)*kz); % column 1:no
	te1 = repmat(te1,2,1);
	tc = 1 ./ (1 - SM1(:,2,2) .* SM2(:,1,1) .* (te1 .* te1));
	Va(1:2*no,1) = ( V(:,1).*SM1(:,2,1) + V(:,2).*SM2(:,1,2).*SM1(:,2,2).*te1 ) .* tc; % e,h+
	Va(2*no+1:4*no,1) = ( V(:,1).*SM1(:,2,1).*SM2(:,1,1).*te1 + V(:,2).*SM2(:,1,2) ) .* tc; % e,h-
	Va = repmat(Va(:,1),1,ns);

	tM = exp(1i*kz.*(zs+0.5*kh)); % matrix no,ns
	Va(1:no,:) = Va(1:no,:) .* tM; % e+
	Va(no+1:2*no,:) = Va(no+1:2*no,:) .* tM; % h+
	tM = exp(1i*kz.*(0.5*kh-zs)); % matrix no,ns
	Va(2*no+1:3*no,:) = Va(2*no+1:3*no,:) .* tM;
	Va(3*no+1:4*no,:) = Va(3*no+1:4*no,:) .* tM;

	U = zeros(3*no, ns);
	
	tM = (Va(1:no,:) + Va(2*no+1:3*no,:)) ./ kxy;
	U(1:no,:) = tM .* ky;
	U(no+1:2*no,:) = -tM .* kx;
	tM = (Va(3*no+1:4*no,:) - Va(no+1:2*no,:)) .* ((1/epsb) * kz ./ kxy);
	U(1:no,:) = U(1:no,:) + tM .* kx;
	U(no+1:2*no,:) = U(no+1:2*no,:) + tM .* ky;
	U(2*no+1:3*no,:) = (Va(3*no+1:4*no,:) + Va(no+1:2*no,:)) .* ((1/epsb) * kxy);
	
	U = reshape(U, [3*no*ns,1]);
end

%%% END OF FILE %%%