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
%% calculate vector of diffraction efficiencies
function [Veff] = gsm_efficiency(Vin, Vout, method, kz1, kz2, eps1, eps2)
	no = method.no(1)*method.no(2);

	ib1 = 1:no;
	ib2 = no+1:2*no;
	
	Veff = zeros(2*no,2);
		% accumulate incident and diffracted field power for each diffraction order
	Pin = sum( abs(Vin(ib1,1).*Vin(ib1,1)).*real(kz1) + abs(Vin(ib1,2).*Vin(ib1,2)).*real(kz2) ) ...
			+ sum( abs(Vin(ib2,1).*Vin(ib2,1)).*real(kz1/eps1) + abs(Vin(ib2,2).*Vin(ib2,2)).*real(kz2/eps2) );
	Veff(ib1,1) = abs(Vout(ib1,1).*Vout(ib1,1)).*real(kz1);
	Veff(ib1,2) = abs(Vout(ib1,2).*Vout(ib1,2)).*real(kz2);
	Veff(ib2,1) = abs(Vout(ib2,1).*Vout(ib2,1)).*real(kz1/eps1);
	Veff(ib2,2) = abs(Vout(ib2,2).*Vout(ib2,2)).*real(kz2/eps2);

	if abs(Pin) > 1e-15
		Veff = (1/Pin) * Veff;
	else
		Veff = 0.5*Veff;
	end
end


