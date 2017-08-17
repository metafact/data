i = find(A.Time == 30/365)

Unik = unique(A.Settlement(i))
for n = 1:length(Unik)

	%PREPARING MATRIX FOR EACH UNIQUE SETTLEMENT TIME
	clear B C D E F
	B = A(i,:);
	j = find(contains(B.Settlement, Unik(n)))

	C = B(j,:);
	C.T = C.Time;
	C.Time = [];
	C.Maturity = [];
	C.Properties.VariableNames{'Vols'} = 'sigmaCall'
	C.Properties.VariableNames{'F'} = 'S'
	C.Properties.VariableNames{'CallPrice'} = 'C'
	E = [table(C.K) C];
	E.K = []
	E.Properties.VariableNames{'Var1'} = 'K'
	E.rf = repmat(0.05, length(E.K), 1);
	F = sortrows(E);
	D = F;


	%START TO INTERPOLATE IMPLIED VOLATILIES clear approxCallPDFs fineK pdfK newC sigmaSplSm

	T0 = unique(D.T)
	extrapThresh = 0.01;
	fineK = linspace(min(D.K)-extrapThresh, ...
					 max(D.K)+extrapThresh, 500).'

	sigmaSplSm = csaps(D.K, D.sigmaCall, 0, fineK)

	newC = NaN(size(sigmaSplSm))

	%ESTIMATE CALL PRICES FROM INTERPOLATED IVs
	S = D.S(1);
	rf = D.rf(1);
	for k = 1:numel(T0)
		newC(:, k) = blsprice(S, fineK, rf, T0(k), sigmaSplSm(:, k));

	end


	%% Estimate the implied densities by approximating derivatives.
	% Approximate the first derivative, at each distinct expiry time.
	% We use the discrete approximation to the first derivative.
	dK = diff(fineK);
	Cdash = diff(newC) ./ repmat(dK, 1, size(newC, 2))

	% Approximate the second derivatives.
	d2K = dK(2:end)
	Cddash = diff(Cdash) ./ repmat(d2K, 1, size(Cdash, 2))

	%% Use the approximate derivatives to estimate each density function.
	% We remove any negative values prior to the estimation process.
	% Preallocate space for the PDF approximations.
	approxCallPDFs = NaN(size(Cddash));

	% Set any negative values to zero.
	Cddash(Cddash < 0) = 0;

	for k = 1:size(Cddash, 2)
		approxCallPDFs(:, k) = exp(rf * T0(k)) * Cddash(:, k);
	end % for

	%% Plot the functions approximated using this technique.
	pdfK = fineK(3:end);
	%plot(pdfK, approxCallPDFs(:, k), 'LineWidth', 2)


	%--------------------%
	%GEV


	%% Fit the right tail with GEV distribution.
	 try
        i = 2
        p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
        p = p * exp(D.rf(1)*T0) +1
        while p <0.84
            i = i+1
            p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
            p = p * exp(D.rf(1)*T0) +1
        end
	end

	if i > 498
		%continue

	else
		a0R = p;
		a1R = 0.995
		Ka0R = fineK(i);
		Ka1R = pdfK(end);
		fKa0R = approxCallPDFs(i-2)
		fKa1R = approxCallPDFs(end)

		%FIND THREE GEV PARAMETERS
		start = [1000 100 -2];
		options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e6, 'TolX', 1e-10);
		[betaR feval] = fminsearch(@(b) FitGEVRightTail(b,a0R,a1R,Ka0R,Ka1R,fKa1R,fKa0R), start, options);
		Mean  = betaR(1);
		sigma = betaR(2);
		phi   = betaR(3);

		% Check how close the optimization conditions are.
		checkCDF_a0 = CDF_GEV(Ka0R,Mean,sigma,phi) - a0R;
		checkPDF_a0 = PDF_GEV(Ka0R,Mean,sigma,phi) - fKa0R;
		checkPDF_a1 = PDF_GEV(Ka1R,Mean,sigma,phi) - fKa1R;
		Error = [checkCDF_a0 checkPDF_a0 checkPDF_a1];

		%ESTIMATE GEV PDFs
		clear GEV_right
		dK = fineK(2)-fineK(1)
		K3 = [pdfK(end):dK:max(D.K)+30];
		dK = fineK(2)-fineK(1)
		GEV_right = gevpdf(K3,phi,sigma,Mean)
		nyK = [pdfK', K3(2:end)]
		nyRND = [approxCallPDFs', GEV_right(2:end)]
	end

	%% Fit the left tail with GEV distribution.

	try
		i = 2
		p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
		p = p * exp(D.rf(1)*T0) +1
		while p <0.05  % Mojet ot desyati procentov nachatsya DAJE??
			i = i+1
			p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
			p = p * exp(D.rf(1)*T0) +1
		end

	end

	if i < 5
		%continue
	else
		a0L = p
		a1L = 0
		Ka0L = fineK(i)
		fKa0L = approxCallPDFs(i-2)

		Ka1L = fineK(3)
		fKa1L = approxCallPDFs(1)
		%%--------------------

		% Find the GEV parameters for the left tail
		% start = betaR;
		start = [1000 100 -2];
		options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e6, 'TolX', 1e-10, 'TolFun', 1e-10);
		[betaL feval] = fminsearch(@(b) FitGEVLeftTail(b,a0L,a1L,Ka0L,Ka1L,fKa1L,fKa0L), start, options);
		Mean  = betaL(1);
		sigma = betaL(2);
		phi   = betaL(3);

		checkCDF_a0 = CDF_GEV(-Ka0L,Mean,sigma,phi) - (1-a0L);
		checkPDF_a0 = PDF_GEV(-Ka0L,Mean,sigma,phi) - fKa0L;
		checkPDF_a1 = PDF_GEV(-Ka1L,Mean,sigma,phi) - fKa1L;
		Error = [checkCDF_a0 checkPDF_a0 checkPDF_a1];

		clear GEV_left
		dK = fineK(2)-fineK(1)
		K4 = linspace(12, pdfK(1), 500)
		GEV_left = gevpdf(-K4,phi,sigma,Mean)
		nyK = [K4, nyK(2:end)]
		nyRND = [GEV_left, nyRND(2:end)]

	end


end
