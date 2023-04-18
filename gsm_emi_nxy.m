%{
Copyright © 2023 Alexey A. Shcherbakov. All rights reserved.
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
function [U, W] = gsm_emi_nxy(V, method)
	no = method.no; % orders
	ns = method.nsl; % slices
	kh = method.kh; % grating depth
	CM_eps = method.CM_eps;
	CM_nxy = method.CM_nxy;
	epsb = method.eps_b;
	zweight = method.zweight;

	NO = no(1)*no(2);
	
	V = reshape(V, [3*NO, ns]);
	U = zeros(3*NO, ns);
	W = zeros(3*NO, ns);

	for is = 1:ns
		U1 = v_expand(V(1:NO, is), no(1), no(2));
		U1 = fft2(U1);
		U2 = ifft2(CM_nxy{1} .* U1); % [nx2]*V
		U3 = ifft2(CM_nxy{3} .* U1); % [nxy]*V
		U1 = ifft2(CM_eps{2} .* U1); % [1/e]*V
		V1 = v_contract(U1, no(1), no(2));
		V2 = v_contract(U2, no(1), no(2));
		V3 = v_contract(U3, no(1), no(2));
		U(1:NO, is) = V2 - V1;
		U(NO+1:2*NO, is) = V3;
		W(1:NO, is) = V1;
		U1 = v_expand(V1, no(1), no(2));
		U1 = ifft2(CM_eps{1} .* fft2(U1)); % [e]*U1
		V1 = v_contract(U1, no(1), no(2));
		U(1:NO, is) = U(1:NO, is) + V1;
		U1 = v_expand(V1, no(1), no(2));
		U1 = fft2(U1);
		U2 = ifft2(CM_nxy{3} .* U1); % [nxy]*U1
		U1 = ifft2(CM_nxy{1} .* U1); % [nx2]*U1
		U(1:NO, is) = U(1:NO, is) - v_contract(U1, no(1), no(2));
		U(NO+1:2*NO, is) = U(NO+1:2*NO, is) - v_contract(U2, no(1), no(2));

		U1 = v_expand(V(NO+1:2*NO, is), no(1), no(2));
		U1 = fft2(U1);
		U2 = ifft2(CM_nxy{2} .* U1); % [ny2]*V
		U3 = ifft2(CM_nxy{3} .* U1); % [nxy]*V
		U1 = ifft2(CM_eps{2} .* U1); % [1/e]*V
		V1 = v_contract(U1, no(1), no(2));
		V2 = v_contract(U2, no(1), no(2));
		V3 = v_contract(U3, no(1), no(2));
		U(NO+1:2*NO, is) = U(NO+1:2*NO, is) + V2 - V1;
		U(1:NO, is) = U(1:NO, is) + V3;
		W(NO+1:2*NO, is) = V1;
		U1 = v_expand(V1, no(1), no(2));
		U1 = ifft2(CM_eps{1} .* fft2(U1)); % [e]*U1
		V1 = v_contract(U1, no(1), no(2));
		U(NO+1:2*NO, is) = U(NO+1:2*NO, is) + V1;
		U1 = v_expand(V1, no(1), no(2));
		U1 = fft2(U1);
		U2 = ifft2(CM_nxy{3} .* U1); % [nxy]*U1
		U1 = ifft2(CM_nxy{2} .* U1); % [ny2]*U1
		U(NO+1:2*NO, is) = U(NO+1:2*NO, is) - v_contract(U1, no(1), no(2));
		U(1:NO, is) = U(1:NO, is) - v_contract(U2, no(1), no(2));

		U1 = v_expand(V(2*NO+1:3*NO, is), no(1), no(2));
		U1 = ifft2(CM_eps{1} .* fft2(U1));
		V1 = v_contract(U1, no(1), no(2));
		W(2*NO+1:3*NO, is) = V1;
		U(2*NO+1:3*NO, is) = V1 - V(2*NO+1:3*NO, is);
		
		U(:, is) = (-0.5*1i*epsb*zweight) * U(:, is);
	end
	
	W = reshape(W, [3*NO*ns,1]);
end

%%% END OF FILE %%%