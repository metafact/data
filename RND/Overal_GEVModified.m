snitt = []
avvik = []
skew = []
kurtos = []
ks = {}
RNDs = {}
nyK = []
nyRND = []

%nanEr = []

m = find(A.Time == 30/365)

Unik = unique(A.Settlement(m))
for n = 1:length(Unik)
    
    %PREPARING MATRIX FOR EACH UNIQUE SETTLEMENT TIME
    clear B C D E F
    B = A(m,:);
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
    
    %SELECT ONLY AT and OUT-OF_THE MONEY
    ut = find(D.K>=D.S);
    D = D(ut,:);
    
    
    %START TO INTERPOLATE IMPLIED VOLATILIES clear approxCallPDFs fineK pdfK newC sigmaSplSm
    clear T0 fineK approxCallPDFs pdfK sigmaSplSm newC Cdash Cddash d2K
    T0 = unique(D.T)
    rf = D.rf(1)
    extrapThresh = 0.01;
    fineK = linspace(min(D.K)-extrapThresh, ...
        max(D.K)+extrapThresh, 500).'
    try
%                     sigmaSplSm = csaps(D.K, D.sigmaCall, 0, fineK)
                pol = polyfit(D.K, D.sigmaCall,4)
                sigmaSplSm = polyval(pol, fineK)
                plot(fineK, sigmaSplSm)
    end
    if ~exist("sigmaSplSm")
        snitt = [snitt NaN]
        avvik = [avvik NaN]
        skew = [skew NaN]
        kurtos = [kurtos NaN]
        continue
    end
    
%          sp = spap2(1,4,D.K, D.sigmaCall)
%          sigmaSplSm = fnval(sp, fineK)
%          fineK = fineK(2:end-1)
%          sigmaSplSm = sigmaSplSm(2:end-1)
    
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
    %start alpha=0.92
    [indsR, probsR] = rightFit(fineK, newC, rf, T0, 0.92, 1)
    
    if ~isnan(indsR) & ~isnan(probsR)
        a0R = probsR(1);
        a1R = 0.995
        Ka0R = fineK(indsR(1));
        Ka1R = pdfK(indsR(2)-2);
        fKa0R = approxCallPDFs(indsR(1)-2)
        fKa1R = approxCallPDFs(indsR(2)-2)
        
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
        clear GEV_right K3 nyK nyRND dK
        dK = fineK(2)-fineK(1)
        K3 = [pdfK(end):dK:max(D.K)+30];
       GEV_right = gevpdf(K3,phi,sigma,Mean)
        
        % 		nyK = [pdfK', K3(2:end)]
        % 		nyRND = [approxCallPDFs', GEV_right(2:end)]
        
        %TEST IF GEV_RIGHT WAS CORRECTLY ESTIMATED, IE. IS BIGGER THAN THE
        %MOST RIGHT POINT OF ESTIMATED RND
        zr = find(GEV_right(2:end)>approxCallPDFs(end))
        if isempty(zr) & GEV_right(1) ~= 0 & max(GEV_right) <= 1 & GEV_right(2) ~= 0
            %if GEV_right(1) ~= 0 & max(GEV_right) <= 1 & GEV_right(2) ~= 0
            nyK = [pdfK(1:end-1)', K3]
            nyRND = [approxCallPDFs(1:end-1)', GEV_right ]
            % 				ks{n} = nyK
            % 				RNDs{n} = nyRND
            
        else
            ks{n} = []
            RNDs{n} = []
            
            %nanEr = [nanEr n]
            snitt = [snitt NaN]
            avvik = [avvik NaN]
            skew = [skew NaN]
            kurtos = [kurtos NaN]
            continue
        end
    else
        ks{n} = []
        RNDs{n} = []
        
        
        %             nanEr = [nanEr n]
        snitt = [snitt NaN]
        avvik = [avvik NaN]
        skew = [skew NaN]
        kurtos = [kurtos NaN]
        
        continue
    end
    
    
    
    %% Fit the left tail with GEV distribution.
    
    
    
    
    [inds, probs] = leftFit(fineK, newC, rf, T0, 0.05, 1)
    
    if ~isnan(inds) & ~isnan(probs)
        a0L = probs(1)
        a1L = 0
        Ka0L = fineK(inds(1))
        fKa0L = approxCallPDFs(inds(1)-2)
        
        Ka1L = fineK(inds(2))
        fKa1L = approxCallPDFs(inds(2)-2)
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
        
        clear GEV_left K4 dK
        dK = fineK(2)-fineK(1)
        K4 = linspace(5, pdfK(1), 500)
        GEV_left = gevpdf(-K4,phi,sigma,Mean)
        
        %TEST IF GEV_LEFT WAS CORRECTLY ESTIMATED, IE. IS BIGGER THAN THE
        %MOST LEFT POINT OF ESTIMATED RND
        zl = find(GEV_left(1:end-2)>approxCallPDFs(1))
        if isempty(zl) & GEV_left(end) ~= 0 & max(GEV_left) <= 1 & GEV_left(end-1) ~= 0
            %if  GEV_left(end) ~= 0 & max(GEV_left) <= 1 & GEV_left(end-1) ~= 0
            %
            % 				nyK = [K4, pdfK(2:end)']
            % 				nyRND = [GEV_left, approxCallPDFs(2:end)']
            nyK = [K4, nyK(2:end)]
            nyRND = [GEV_left, nyRND(2:end)]
            % 				ks{n} = nyK
            % 				RNDs{n} = nyRND
            
        else
            ks{n} = []
            RNDs{n} = []
            
            %                 nanEr = [nanEr n]
            snitt = [snitt NaN]
            avvik = [avvik NaN]
            skew = [skew NaN]
            kurtos = [kurtos NaN]
            continue
        end
    else
        ks{n} = []
        RNDs{n} = []
        
        %             nanEr = [nanEr n]
        snitt = [snitt NaN]
        avvik = [avvik NaN]
        skew = [skew NaN]
        kurtos = [kurtos NaN]
        continue
    end
    ks{n} = nyK
    RNDs{n} = nyRND
    
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
end


%%%%%NIJE TESTOVIY VARIANT

% M3 = (nyK -avrg).^3
% M3 = trapz(nyK, (M3.*nyRND) )
%
% skew = [skew M3]
%
% M4 = (nyK -avrg).^4
% M4 = trapz(nyK, (M4.*nyRND) )
%
% kurtos = [kurtos M4]





%Initial value t=1 and intially value alpha must be 0.05


% xticklabels(dater)
% tel = find(~isnan(skew))
%
%
% %%%%%%%%%%%%%%%%%%%%
% %LOAD FUT.CSV FILE
% fid = fopen('fut.csv')
%       fut = textscan(fid,'%s %s %f %f %f %f %f %f', 'Delimiter',',')
% fclose(fid)
%
%
%  %PREPARE DATES FORMAT
% Mat = regexprep(fut{1}(1:end-1), 'Z', '/11/19');
% Mat1 = regexprep(Mat, 'K', '/04/19');
% Mat2 = regexprep(Mat1, 'G', '/01/19');
% Mat3 = regexprep(Mat2, 'H', '/02/19');
% Mat4 = regexprep(Mat3, 'J', '/03/19');
% Mat5 = regexprep(Mat4, 'M', '/05/19');
% Mat6 = regexprep(Mat5, 'N', '/06/19');
% Mat7 = regexprep(Mat6, 'Q', '/07/19');
% Mat8 = regexprep(Mat7, 'U', '/08/19');
% Mat9 = regexprep(Mat8, 'V', '/09/19');
% Mat10 = regexprep(Mat9, 'X', '/10/19');
% Mat11 = regexprep(Mat10, 'F', '/12/19');
%
% Mat11 = char(Mat11)
% futpr = Mat11(:,3:end)
%
%
% %PREPARING CELL ARRAY FOR ESTIMATION OF EMPIRICAL MOMENTS
% clear fc
% MatNum = datenum(futpr)
% Unik2 = Unik(tel)
% for ii = 1:length(Unik2)
%     DatNum = datenum(Unik2(ii,:))
%     DatNum1M = DatNum + 30
%     g = find(MatNum == DatNum1M)
%     if ~isempty(g)
%         fc{ii} = fut{6}(g)
%     else
%         fc{ii} = NaN
%     end
%
% end
%
% %ESTIMATE EMPIRICAL MOMENTS
% clear snittE avvikE skewE kurtosE std2 std
% snittE = []
% avvikE = []
% skewE = []
% kurtosE = []
%
% for jj = 1:length(fc)
%     mn = mean(fc{jj})
%     snittE = [snittE, mn]
%
%     av = std2(fc{jj})
%     avvikE = [avvikE, av]
%
%     sk = skewness(fc{jj})
%     skewE = [skewE, sk]
%
%     kr = kurtosis(fc{jj})
%     kurtosE = [kurtosE, kr]
% end
% moments = [datenum(Unik(tel)),snitt(tel)', avvik(tel)', skew(tel)', kurtos(tel)', ...
% snittE', avvikE', skewE', kurtosE']
%
% nmom = sortrows(moments,1)
%
%  plot(nmom(:,1), (nmom(:,3)./T0).^0.5,nmom(:,1), nmom(:,7))
%  plot(nmom(:,1), nmom(:,2),nmom(:,1), nmom(:,6))
%
%  plot(nmom(:,1), nmom(:,4),nmom(:,1), nmom(:,8))
%  plot(nmom(:,1), nmom(:,5),nmom(:,1), nmom(:,9))
%
% plot(nmom(:,1), nmom(:,5))
% xtickangle(90)
% xticks([1:4:68])
% datetick('x',20)
% grid on
% datetick('x',20)
% legend('Risk Neutral', 'Empirical')
% title('Skewness, 1 M')
% ylabel('USD')


% %TESTING PLOTS
ze = find(~cellfun(@isempty, ks))
%for i = 46:round(length(ze)/2)
for i = 1:length(ze)
    figure(i)
    plot(ks{ ze(i) }, RNDs{ ze(i) })
end

ze = find(~isnan(avvik))
length(ze)
M5 = table(Unik, snitt', avvik', skew', kurtos' )
save('M5.mat', 'M5')

% for i = 40:53
%     close(figure(i))
%
% end