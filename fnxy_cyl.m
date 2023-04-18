%% implementation

function FN = fnxy_cyl(xno, yno, kgx, kgy, kr1, kr0, kr2, cx, cy)
	FN = cellmat(3,1,2*yno-1,2*xno-1); % cos^2, sin^2, cos*sin

	nc = numel(kr0);
	if nc ~= numel(kr1) || nc ~= numel(kr2) || nc ~= numel(cx) || nc ~= numel(cy)
		error('incorrect parameters in function fnxy_cyl');
	end

	ix = linspace(-xno+1,xno-1,2*xno-1);
	iy = linspace(-yno+1,yno-1,2*yno-1);
	[IX,IY] = meshgrid(ix,iy);
	
	for i = 1:nc
			% integration parameters:
		CC = [(kr1(i)+kr2(i)+kr0(i))*(kr2(i)-kr1(i))/12, ...
					kr0(i)*(kr2(i)-kr1(i))/((kr0(i)-kr1(i))*(kr2(i)-kr0(i))), -kr1(i)/(kr0(i)-kr1(i)), -kr2(i)/(kr2(i)-kr0(i)), ...
					-2*(kr2(i)-kr1(i))/((kr0(i)-kr1(i))*(kr2(i)-kr0(i))), 2/(kr0(i)-kr1(i)), 2/(kr2(i)-kr0(i))];
				
		M = cellmat(3,1,2*yno-1,2*xno-1);

		for m = 0:xno-1
			for n = -yno+1:yno-1
				fun_c2 = @(x) fun_int_c2(x, kr0(i), kr1(i), kr2(i), CC, m*kgx, n*kgy);
				fun_s2 = @(x) fun_int_s2(x, kr0(i), kr1(i), kr2(i), CC, m*kgx, n*kgy);
				fun_sc = @(x) fun_int_sc(x, kr0(i), kr1(i), kr2(i), CC, m*kgx, n*kgy);
				M{1}(yno+n,xno+m) = quadgk(fun_c2, 0, pi, 'RelTol', 1e-12);%, 'MaxIntervalCount', 1e4*(1+abs(m)+abs(n)));
				M{2}(yno+n,xno+m) = quadgk(fun_s2, 0, pi, 'RelTol', 1e-12);%, 'MaxIntervalCount', 1e4*(1+abs(m)+abs(n)));
				M{3}(yno+n,xno+m) = quadgk(fun_sc, 0, pi, 'RelTol', 1e-12);%, 'MaxIntervalCount', 1e4*(1+abs(m)+abs(n)));
			end
		end

		EXP = exp(-(1i*2*pi*cx(i))*IX - (1i*2*pi*cy(i))*IY);

		for i = 1:3
			M{i}(yno:-1:1, xno-1:-1:1) = M{i}(yno:2*yno-1,xno+1:2*xno-1);
			M{i}(yno+1:2*yno-1, xno-1:-1:1) = M{i}(yno-1:-1:1,xno+1:2*xno-1);
			FN{i} = FN{i} + (0.5*kgx*kgy/pi/pi) * M{i} .* EXP;
		end
	end % end of loop over cylinders
end

function f = fun_int_c2(phi, kr0, kr1, kr2, CC, mgx, ngy)
	kappa = mgx*cos(phi) + ngy*sin(phi);
	ikappa = 1 ./ kappa;
	cosr0 = cos(kr0 * kappa);
	cosr1 = cos(kr1 * kappa);
	cosr2 = cos(kr2 * kappa);
	sinr0 = sin(kr0 * kappa);
	sinr1 = sin(kr1 * kappa);
	sinr2 = sin(kr2 * kappa);

	ind0 = abs(kappa) < 1e-7;
	ind1 = ~ind0;
	f = 0*kappa;

	if sum(ind1,'all') > 0
		f(ind1) = ( (ikappa(ind1).^2) .* ( CC(2) * cosr0(ind1) + CC(3) * cosr1(ind1) + CC(4) * cosr2(ind1) ) ...
							+ (ikappa(ind1).^3) .* ( CC(5) * sinr0(ind1) + CC(6) * sinr1(ind1) + CC(7) * sinr2(ind1) ) ) ...
							.* ((cos(phi(ind1))).^2);
	end
	if sum(ind0,'all') > 0
		f(ind0) = f(ind0) + CC(1);
	end
end

function f = fun_int_s2(phi, kr0, kr1, kr2, CC, mgx, ngy)
	kappa = mgx*cos(phi) + ngy*sin(phi);
	ikappa = 1 ./ kappa;
	cosr0 = cos(kr0 * kappa);
	cosr1 = cos(kr1 * kappa);
	cosr2 = cos(kr2 * kappa);
	sinr0 = sin(kr0 * kappa);
	sinr1 = sin(kr1 * kappa);
	sinr2 = sin(kr2 * kappa);

	ind0 = abs(kappa) < 1e-7;
	ind1 = ~ind0;
	f = 0*kappa;

	if sum(ind1,'all') > 0
		f(ind1) = ( (ikappa(ind1).^2) .* ( CC(2) * cosr0(ind1) + CC(3) * cosr1(ind1) + CC(4) * cosr2(ind1) ) ...
							+ (ikappa(ind1).^3) .* ( CC(5) * sinr0(ind1) + CC(6) * sinr1(ind1) + CC(7) * sinr2(ind1) ) ) ...
							.* ((sin(phi(ind1))).^2);
	end
	if sum(ind0,'all') > 0
		f(ind0) = f(ind0) + CC(1);
	end
end

function f = fun_int_sc(phi, kr0, kr1, kr2, CC, mgx, ngy)
	kappa = mgx*cos(phi) + ngy*sin(phi);
	ikappa = 1 ./ kappa;
	cosr0 = cos(kr0 * kappa);
	cosr1 = cos(kr1 * kappa);
	cosr2 = cos(kr2 * kappa);
	sinr0 = sin(kr0 * kappa);
	sinr1 = sin(kr1 * kappa);
	sinr2 = sin(kr2 * kappa);

	ind0 = abs(kappa) < 1e-7;
	ind1 = ~ind0;
	f = 0*kappa;

	if sum(ind1,'all') > 0
		f(ind1) = ( (ikappa(ind1).^2) .* ( CC(2) * cosr0(ind1) + CC(3) * cosr1(ind1) + CC(4) * cosr2(ind1) ) ...
							+ (ikappa(ind1).^3) .* ( CC(5) * sinr0(ind1) + CC(6) * sinr1(ind1) + CC(7) * sinr2(ind1) ) ) ...
							.* ((sin(phi(ind1))).*(cos(phi(ind1))));
	end
end

%%% end of file %%%