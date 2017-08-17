snitt = []
avvik = []
skew = []
kurtos = []
%%%%%%%%

teller = 0;


m = find(Uniq.MatRatio == 30/365)

Unik = unique(Uniq.Date(m))
for n = 1:100%length(Unik)

	%PREPARING MATRIX FOR EACH UNIQUE SETTLEMENT TIME
	clear B C D E F
	B = Uniq(m,:); %Table with only 30 days to maturity
	j = find(contains(B.Date, Unik(n)))

	C = B(j,:);
	C.T = C.MatRatio;
	C.MatRatio = [];
	C.Maturity = [];
	C.Properties.VariableNames{'Vols'} = 'sigmaCall'
	C.Properties.VariableNames{'S0'} = 'S'
	C.Properties.VariableNames{'OptPrice'} = 'C'
	E = [table(C.K) C];
	E.K = []%l.ll ???
	E.Properties.VariableNames{'Var1'} = 'K'
	E.rf = E.Rate;
	F = sortrows(E);
	D = F;

    
   %%%%%%%%%%%%%%%
   if max(D.sigmaCall) > 10
       teller = teller + 1;
       continue
   end

	%START TO INTERPOLATE IMPLIED VOLATILIES clear approxCallPDFs fineK pdfK newC sigmaSplSm
    clear T0 fineK approxCallPDFs pdfK sigmaSplSm newC Cdash Cddash d2K
	T0 = unique(D.T)
	extrapThresh = 0.01;
	fineK = linspace(min(D.K)-extrapThresh, ...
					 max(D.K)+extrapThresh, 500).'

	sigmaSplSm = csaps(D.K, D.sigmaCall, 0.001, fineK)
    %plot(D.K, D.sigmaCall, '-ro', fineK, sigmaSplSm, '-x')

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
    figure(n)
	pdfK = fineK(3:end);
	%plot(pdfK, approxCallPDFs(:, k), 'LineWidth', 2)

end
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
    
	if i > 500
        snitt = [snitt NaN]
        avvik = [avvik NaN]
        skew = [skew NaN]
        kurtos = [kurtos NaN]
		continue

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
        snitt = [snitt NaN]
        avvik = [avvik NaN]
        skew = [skew NaN]
        kurtos = [kurtos NaN]
		continue
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
    avrg = trapz( nyK, (nyK.*nyRND) ) 
    snitt = [snitt avrg]
    
    std = (nyK -avrg).^2
    std2 = trapz( nyK, (std.*nyRND) )
    avvik = [avvik std2]
    
    M2 = std2
    
    M3 = (nyK -avrg).^3
    M3 = trapz(nyK, (M3.*nyRND) ) 
    skjevhet = M3/(M2^(3/2))
    skew = [skew skjevhet]
    
    M4 = (nyK -avrg).^4
    M4 = trapz(nyK, (M4.*nyRND) )
    kurts = M4/(M2^(2))
    kurtos = [kurtos kurts]

% M3 = (nyK -avrg).^3
%     M3 = trapz(nyK, (M3.*nyRND) ) 
%    
%     skew = [skew M3]
%     
%     M4 = (nyK -avrg).^4
%     M4 = trapz(nyK, (M4.*nyRND) )
%     
%     kurtos = [kurtos M4]



end

% xticklabels(dater)
tel = find(~isnan(skew))


%%%%%%%%%%%%%%%%%%%%
%LOAD FUT.CSV FILE
fid = fopen('fut.csv')
      fut = textscan(fid,'%s %s %f %f %f %f %f %f', 'Delimiter',',')
fclose(fid)
 
 
 %PREPARE DATES FORMAT
Mat = regexprep(fut{1}(1:end-1), 'Z', '/11/19');
Mat1 = regexprep(Mat, 'K', '/04/19');
Mat2 = regexprep(Mat1, 'G', '/01/19');
Mat3 = regexprep(Mat2, 'H', '/02/19');
Mat4 = regexprep(Mat3, 'J', '/03/19');
Mat5 = regexprep(Mat4, 'M', '/05/19');
Mat6 = regexprep(Mat5, 'N', '/06/19');
Mat7 = regexprep(Mat6, 'Q', '/07/19');
Mat8 = regexprep(Mat7, 'U', '/08/19');
Mat9 = regexprep(Mat8, 'V', '/09/19');
Mat10 = regexprep(Mat9, 'X', '/10/19');
Mat11 = regexprep(Mat10, 'F', '/12/19');

Mat11 = char(Mat11)
futpr = Mat11(:,3:end)


%PREPARING CELL ARRAY FOR ESTIMATION OF EMPIRICAL MOMENTS
clear fc
MatNum = datenum(futpr)
Unik2 = Unik(tel)
for ii = 1:length(Unik2)
    DatNum = datenum(Unik2(ii,:))
    DatNum1M = DatNum + 30
    g = find(MatNum == DatNum1M)
    if ~isempty(g)
        fc{ii} = fut{6}(g)
    else    
        fc{ii} = NaN
    end 
    
end 

%ESTIMATE EMPIRICAL MOMENTS
clear snittE avvikE skewE kurtosE std2 std
snittE = []
avvikE = []
skewE = []
kurtosE = []

for jj = 1:length(fc)
    mn = mean(fc{jj})
    snittE = [snittE, mn]
    
    av = std2(fc{jj})
    avvikE = [avvikE, av]
    
    sk = skewness(fc{jj})
    skewE = [skewE, sk]
    
    kr = kurtosis(fc{jj})
    kurtosE = [kurtosE, kr]
end
moments = [datenum(Unik(tel)),snitt(tel)', avvik(tel)', skew(tel)', kurtos(tel)', ...
snittE', avvikE', skewE', kurtosE']
 
nmom = sortrows(moments,1)

 plot(nmom(:,1), (nmom(:,3)./T0).^0.5,nmom(:,1), nmom(:,7))
 plot(nmom(:,1), nmom(:,2),nmom(:,1), nmom(:,6))
 
 plot(nmom(:,1), nmom(:,4),nmom(:,1), nmom(:,8))
 plot(nmom(:,1), nmom(:,5),nmom(:,1), nmom(:,9))
 
plot(nmom(:,1), nmom(:,5))
xtickangle(90)
xticks([1:4:68])
datetick('x',20)
grid on
datetick('x',20)
legend('Risk Neutral', 'Empirical')
title('Skewness, 1 M')
ylabel('USD')