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
%% add incident field amplitudes to the calculated diffracted field amplitudes:
function U = gsm_add_inc(Vin, Vout, method)
	no = method.no(1)*method.no(2); % orders
	kh = method.kh;
	kz = method.kz;
	SM1 = method.SMsubb;
	SM2 = method.SMbsup;
	
	SM = mul_SM(SM1, smatrix_layer(kz, kh, 'TEM'));
	SM = mul_SM(SM, SM2);
	
	U = Vout;
	U(:,1) = U(:,1) + SM(:,1,1).*Vin(:,1) + SM(:,1,2).*Vin(:,2);
	U(:,2) = U(:,2) + SM(:,2,1).*Vin(:,1) + SM(:,2,2).*Vin(:,2);
end

%%% END OF FILE %%%