function [avrg, std2, skjevhet, kurts]= RiNeDe(matrise)
%function stor = RiNeDe(matrise)
clear B C D E F
C = matrise;
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


stor = find(D.sigmaCall < 0 | D.sigmaCall > 5)
D.sigmaCall(stor) = NaN
D.K(stor) = NaN  %DOUBLE CHECK THIS

%%%%%%%%%%%%%%
% if max(D.sigmaCall) > 10 | D.sigmaCall < 0
%     %         Ks(end + 1,:) = {NaN}
%     %         RND(end + 1,:) = {NaN}
%     %
%     %         snitt = [snitt NaN]
%     %         avvik = [avvik NaN]
%     %         skew = [skew NaN]
%     %         kurtos = [kurtos NaN]
%     %
%     %         teller = teller + 1;
%     return
% end

%START TO INTERPOLATE IMPLIED VOLATILIES clear approxCallPDFs fineK pdfK newC sigmaSplSm
clear T0 fineK approxCallPDFs pdfK sigmaSplSm newC Cdash Cddash d2K
T0 = unique(D.T)

%T0 = D.T(1)
extrapThresh = 0.01;
fineK = linspace(min(D.K)-extrapThresh, ...
    max(D.K)+extrapThresh, 500).'

sigmaSplSm = csaps(D.K, D.sigmaCall, 0.007, fineK)

neg = find(sigmaSplSm < 0)
sigmaSplSm(neg) = []
fineK(neg) = []
%     figure(n)
%     plot(D.K, D.sigmaCall, '-ro', fineK, sigmaSplSm, '-x')

newC = NaN(size(sigmaSplSm))

%ESTIMATE CALL PRICES FROM INTERPOLATED IVs
ii = 1; %This is done for choosing NaN interest rate
S = D.S(1);
while isnan(D.rf(ii)) && ii <= numel(D.rf)
    ii = ii + 1
end
rf = D.rf(ii);
if rf < 0
    rf = 0
end

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


%%%

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
% stor = size(pdfK,1)
% figure(n)
% plot(pdfK, approxCallPDFs(:, k), 'LineWidth', 2)
%     end
%--------------------%
%GEV

%% Fit the right tail with GEV distribution.



%Estimate for i = 499
%p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
%if p = 99 or close to it then no need for GEV tail?
i = numel(fineK) - 1 % 499
p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
p = p * exp(D.rf(ii)*T0) +1 %Formula for cumulative distribution %SHOULD NEGATTIVE INTEREST RATES BE SET TO ZERO
if p >= 0.995
    %     tail_skip = tail_skip +1;
    nyK_R = pdfK';
    nyRND_R = approxCallPDFs(:, k)';
    %         figure(n)
    %         plot(pdfK, approxCallPDFs(:, k), 'LineWidth', 2)
    
else
    try
        i = 2
        % "i+1" and "i-1" is taken such that First order
        % derivative/Cumulative probability is exactly on "K(i)"/fineK(i)
        % in the GEV three conditions further down
        p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1)) %First order derivative
        %http://www2.warwick.ac.uk/fac/soc/wbs/subjects/finance/research/wpaperseries/1998/98-95.pdf
        %page 8
        p0 = p * exp(D.rf(ii)*T0) +1
        while p0 <0.87
            i = i+1
            p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
            p0 = p * exp(D.rf(ii)*T0) +1
        end
        j = i;
        p1 = p0;
        while p1 <0.93 && j < numel(fineK) - 2 % 498
            j = j+1
            p = (newC(j+1)-newC(j-1))/ (fineK(j+1)-fineK(j-1))
            p1 = p * exp(D.rf(ii)*T0) +1
        end
    end
    
    if i > numel(fineK) - 2 | i == j
        
        %             Ks(end + 1,:) = {NaN}
        %             RND(end + 1,:) = {NaN}
        avrg = NaN
        std2 = NaN
        skjevhet = NaN
        kurts = NaN
        
        %         snitt = [snitt NaN]
        %         avvik = [avvik NaN]
        %         skew = [skew NaN]
        %         kurtos = [kurtos NaN]
        
        %             tail_skip = tail_skip +1;
        %
        
        return
    else
        
        a0R = p0;
        a1R = 0.995 %Is not used
        Ka0R = fineK(i);
        Ka1R = pdfK(j); %Should it be (i-2)
        %Ka1R = pdfK(end);
        fKa0R = approxCallPDFs(i-2) %try (i-1) here
        fKa1R = approxCallPDFs(j)
        
        
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
        K3 = [pdfK(end):dK:max(D.K)+35];
        dK = fineK(2)-fineK(1)
        GEV_right = gevpdf(K3,phi,sigma,Mean)
        nyK_R = [pdfK', K3(2:end)]
        nyRND_R = [approxCallPDFs', GEV_right(2:end)]
        %end
        %             figure(n)
        %             plot(nyK, nyRND);
    end
end
%end


%% Fit the left tail with GEV distribution.

i = 2
p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
p = p * exp(D.rf(ii)*T0) +1 %Formula for cumulative distribution
if p < 0.005
    
    nyK_L= nyK_R
    nyRND_L = nyRND_R
    
    %     Ks(end + 1,:) = {nyK_L}
    %     RND(end + 1,:) = {nyRND_L}
    %
    [avrg, std2, skjevhet, kurts] = est_moments(nyK_L, nyRND_L,S)
    
    return
    
    %         figure(n)
    %         plot(pdfK, approxCallPDFs(:, k), 'LineWidth', 2)
else
    try
        i = 2
        p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
        p = p * exp(D.rf(ii)*T0) + 1
        while p < 0.09  % Mojet ot desyati procentov nachatsya DAJE??
            i = i+1
            p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
            p = p * exp(D.rf(ii)*T0) +1
        end
        j = 2
        p1 = 0
        while p1 < 0.04 && (i-j) > 0
            j = j + 1
            p1 = (newC(j+1)-newC(j-1))/ (fineK(j+1)-fineK(j-1))
            p1 = p1 * exp(D.rf(ii)*T0) +1
        end
    end
    if i < 5
        %         Ks(end + 1,:) = {NaN}
        %         RND(end + 1,:) = {NaN}
        
        
        %         snitt = [snitt NaN]
        %         avvik = [avvik NaN]
        %         skew = [skew NaN]
        %         kurtos = [kurtos NaN]
        %             figure(n)
        %             plot(nyK, nyRND);
        
        avrg = NaN
        std2 = NaN
        skjevhet = NaN
        kurts = NaN
        
        return
        
    else
        a0L = p
        a1L = 0 %Is not USED????
        Ka0L = fineK(i)
        fKa0L = approxCallPDFs(i-2)
        
        Ka1L = pdfK(j)
        fKa1L = approxCallPDFs(j)
        %%--------------------
        
        % Find the GEV parameters for the left tail
        % start = betaR;
        start = [1000 100 -2];
%         options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e8, 'TolX', 1e-10, 'TolFun', 1e-10);
        options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e6, 'TolX', 1e-10);
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
        %K4 = linspace(10, pdfK(1), 500)
        K4 = [4:dK:pdfK(1)];
        GEV_left = gevpdf(-K4,phi,sigma,Mean)
        nyK_L= [K4, nyK_R(2:end)]
        nyRND_L = [GEV_left, nyRND_R(2:end)]
        
        
        [avrg, std2, skjevhet, kurts] = est_moments(nyK_L, nyRND_L,S)
        
    end
end

    function [avrg, std2, skjevhet, kurts]= est_moments(nyK_L, nyRND_L,S)
        
        avrg = trapz( log(nyK_L/S), (log(nyK_L/S).*nyRND_L))
        
        
        std = (log(nyK_L/S) -avrg).^2
        std2 = trapz( log(nyK_L/S), (std.*nyRND_L) )
        
        %
        M2 = std2
        M3 = (log(nyK_L/S) -avrg).^3
        M3 = trapz(log(nyK_L/S), (M3.*nyRND_L) )
        skjevhet = M3/(M2^(3/2))
        
        %
        M4 = (log(nyK_L/S) -avrg).^4
        M4 = trapz(log(nyK_L/S), (M4.*nyRND_L) )
        kurts = M4/(M2^(2))
        
        %
        % % M3 = (nyK -avrg).^3
        % %     M3 = trapz(nyK, (M3.*nyRND) )
        % %
        % %     skew = [skew M3]
        % %
        % %     M4 = (nyK -avrg).^4
        % %     M4 = trapz(nyK, (M4.*nyRND) )
        % %
        % %     kurtos = [kurtos M4]
        %
        %
    end


end