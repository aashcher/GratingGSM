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
%% calculate amplitude vector of the diffracted field
function [Vout, Veff, balance, rv, U] = gsm_diffraction(incidence, method, grating, bnorm, bini)
	Uinc = gsm_inc(incidence.Vinc, method);
	if ~bnorm
		if ~bini
			[U, ~, ~, ~, rv] = gmres(@(v) multiply(v, method), Uinc, method.restart, method.tol, method.maxit);
		else
			[U, ~, ~, ~, rv] = gmres(@(v) multiply(v, method), Uinc, method.restart, method.tol, method.maxit, [], [], method.Uini);
		end
		Vout = gsm_propb(gsm_emi(U, method), method);
	else
		[U, ~, ~, ~, rv] = gmres(@(v) multiply_nxy(v, method), Uinc, method.restart, method.tol, method.maxit);
		Vout = gsm_propb(gsm_emi_nxy(U, method), method);
	end
	Vout = gsm_add_inc(incidence.Vinc, Vout, method);
	fprintf('gsm gmres: %d iterations, residual = %f\n', numel(rv), rv(end));
	
	kg = [incidence.wavelength/grating.period(1), incidence.wavelength/grating.period(2)];
	[kz1, kz2] = kxyz(method.no, incidence.k_inc, kg, [grating.eps_sub, grating.eps_sup]);
	Veff = gsm_efficiency(incidence.Vinc, Vout, method, kz1, kz2, grating.eps_sub, grating.eps_sup);
	
	balance = abs(1 - sum(real(Veff),'all'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Vout = multiply(Vin, method)
	Vout = Vin - gsm_prop(gsm_emi(Vin, method), method);
end

function Vout = multiply_nxy(Vin, method)
	[U, W] = gsm_emi_nxy(Vin, method);
	Vout = W - gsm_prop(U, method);
end

%%% END OF FILE %%%