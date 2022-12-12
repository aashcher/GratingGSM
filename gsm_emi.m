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
%% calculate emission in slices from the generalized sources:
function W = gsm_emi(V, method)
	no = method.no; % orders
	ns = method.nsl; % slices
	kh = method.kh; % grating depth
	CM_eps = method.CM_eps;
	epsb = method.eps_b;
	zweight = method.zweight;

	NO = no(1)*no(2);
	
	W = reshape(V, [3*NO,ns]);

	for is = 1:ns
		tc = -0.5*1i*epsb*zweight;

		U = v_expand(W(1:NO, is), no(1), no(2));
		U = ifft2(CM_eps{1} .* fft2(U)) - U;
		W(1:NO, is) = tc * v_contract(U, no(1), no(2));

		U = v_expand(W(NO+1:2*NO, is), no(1), no(2));
		U = ifft2(CM_eps{1} .* fft2(U)) - U;
		W(NO+1:2*NO, is) = tc * v_contract(U, no(1), no(2));

		U = v_expand(W(2*NO+1:3*NO, is), no(1), no(2));
		U = U - ifft2(CM_eps{2} .* fft2(U));
		W(2*NO+1:3*NO, is) = tc * v_contract(U, no(1), no(2));
	end
end

%%% END OF FILE %%%