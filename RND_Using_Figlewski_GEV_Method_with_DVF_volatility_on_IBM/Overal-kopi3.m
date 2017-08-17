snitt = []
avvik = []
skew = []
kurtos = []
%%%%%%%%

teller = 0;
tail_skip = 0;
monoton = [];
monoton0 = [];

m = find(Uniq.MatRatio == 30/365)

Unik = unique(Uniq.Date(m))
for n = 1:50%length(Unik)
    
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
    
    sigmaSplSm = csaps(D.K, D.sigmaCall, 0.01, fineK)
%     figure(n)
%     plot(D.K, D.sigmaCall, '-ro', fineK, sigmaSplSm, '-x')
    
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
    
    
%     pos = find(Cdash > 0);
%     disp("Number of positive and nonmonton:")
%     disp(pos)
    % Approximate the second derivatives.
    d2K = dK(2:end)
    Cddash = diff(Cdash) ./ repmat(d2K, 1, size(Cdash, 2))
%     
    for t = 2:length(Cdash)
        if Cdash(t) < Cdash(t-1)
            monoton = [monoton t];
            
        end
    end
    
    %%For not interpolated Vols
%     newC0 = NaN(size(D.sigmaCall))
%     
%     %ESTIMATE CALL PRICES FROM INTERPOLATED IVs
%     for k = 1:numel(T0)
%         newC0(:, k) = blsprice(S, D.K, rf, T0(k), D.sigmaCall(:, k));
%         
%     end
%     dK0 = diff(D.K);
%     
%     Cdash0 = diff(newC0) ./ repmat(dK0, 1, size(newC0, 2))
%     
%     d2K0 = dK0(2:end)
%     Cddash0 = diff(Cdash0) ./ repmat(d2K0, 1, size(Cdash0, 2))
%     
%     for t = 2:length(Cdash0)
%         if Cdash0(t) < Cdash0(t-1)
%             monoton0 = [monoton0 t];
%             
%         end
%     end
    
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
    figure(n)
    plot(pdfK, approxCallPDFs(:, k), 'LineWidth', 2)
    
    %end
    %--------------------%
    %GEV
    
    
    %% Fit the right tail with GEV distribution.
    
    
    
    %Estimate for i = 499
    %p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
    %if p = 99 or close to it then no need for GEV tail?
    i = 499
    p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
    p = p * exp(D.rf(1)*T0) +1
    if p >= 0.99
        %tail_skip = tail_skip +1;
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
            p0 = p * exp(D.rf(1)*T0) +1
            while p0 <0.87
                i = i+1
                p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
                p0 = p * exp(D.rf(1)*T0) +1
            end
            j = i;
            p1 = p0;
            while p1 <0.93 && j < 498
                j = j+1
                p = (newC(j+1)-newC(j-1))/ (fineK(j+1)-fineK(j-1))
                p1 = p * exp(D.rf(1)*T0) +1
            end
        end
        
        if i > 498 | i == j
            snitt = [snitt NaN]
            avvik = [avvik NaN]
            skew = [skew NaN]
            kurtos = [kurtos NaN]
            
            tail_skip = tail_skip +1;
            
            continue
            
        else
            
            a0R = p0;
            a1R = 0.995 %Is not used
            Ka0R = fineK(i);
            Ka1R = pdfK(j);
            %Ka1R = pdfK(end);
            fKa0R = approxCallPDFs(i-2) %try (i-1) here
            fKa1R = approxCallPDFs(j)
            %fKa1R = approxCallPDFs(end)
            
            
            
            
            %     %Estimate for i = 499
            %     %p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
            %     %if p = 99 or close to it then no need for GEV tail?
            %     i = 499
            %     p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
            %     p = p * exp(D.rf(1)*T0) +1
            %     if p >= 0.99
            %         tail_skip = tail_skip +1;
            %     else
            %         try
            %             i = 2
            %             % "i+1" and "i-1" is taken such that First order
            %             % derivative/Cumulative probability is exactly on "K(i)"/fineK(i)
            %             % in the GEV three conditions further down
            %             p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1)) %First order derivative
            %             %http://www2.warwick.ac.uk/fac/soc/wbs/subjects/finance/research/wpaperseries/1998/98-95.pdf
            %             %page 8
            %             p = p * exp(D.rf(1)*T0) +1
            %             while p <0.84
            %                 i = i+1
            %                 p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
            %                 p = p * exp(D.rf(1)*T0) +1
            %             end
            %         end
            %
            %         if i > 500
            %             snitt = [snitt NaN]
            %             avvik = [avvik NaN]
            %             skew = [skew NaN]
            %             kurtos = [kurtos NaN]
            %             continue
            %
            %         else
            %             if i > 497
            %                 i = i-10
            %                 %a0R = p;
            %                 %a1R = 0.995 %Is not used
            %                 Ka0R = fineK(i);
            %                 %Ka1R = pdfK(end);
            %                 fKa0R = approxCallPDFs(i-2) %triy (i-1) here
            %                 %fKa1R = approxCallPDFs(end)
            %                 p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
            %                 p = p * exp(D.rf(1)*T0) +1
            %             else
            %                 Ka0R = fineK(i);
            %                 fKa0R = approxCallPDFs(i-2)
            %             end
            %             a0R = p;
            %             a1R = 0.995 %Is not used
            %             Ka1R = pdfK(end);
            %             fKa1R = approxCallPDFs(end)
            
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
            %end
%             figure(n)
%             plot(nyK, nyRND);
        end
    end
end


%% Fit the left tail with GEV distribution.

% try
%     i = 2
%     p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
%     p = p * exp(D.rf(1)*T0) +1
%     while p <0.05  % Mojet ot desyati procentov nachatsya DAJE??
%         i = i+1
%         p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
%         p = p * exp(D.rf(1)*T0) +1
%     end
%     
% end
% 
% if i < 5
%     snitt = [snitt NaN]
%     avvik = [avvik NaN]
%     skew = [skew NaN]
%     kurtos = [kurtos NaN]
%     continue
% else
%     a0L = p
%     a1L = 0
%     Ka0L = fineK(i)
%     fKa0L = approxCallPDFs(i-2)
%     
%     Ka1L = fineK(3)
%     fKa1L = approxCallPDFs(1)
%     %%--------------------
%     
%     % Find the GEV parameters for the left tail
%     % start = betaR;
%     start = [1000 100 -2];
%     options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e6, 'TolX', 1e-10, 'TolFun', 1e-10);
%     [betaL feval] = fminsearch(@(b) FitGEVLeftTail(b,a0L,a1L,Ka0L,Ka1L,fKa1L,fKa0L), start, options);
%     Mean  = betaL(1);
%     sigma = betaL(2);
%     phi   = betaL(3);
%     
%     checkCDF_a0 = CDF_GEV(-Ka0L,Mean,sigma,phi) - (1-a0L);
%     checkPDF_a0 = PDF_GEV(-Ka0L,Mean,sigma,phi) - fKa0L;
%     checkPDF_a1 = PDF_GEV(-Ka1L,Mean,sigma,phi) - fKa1L;
%     Error = [checkCDF_a0 checkPDF_a0 checkPDF_a1];
%     
%     clear GEV_left
%     dK = fineK(2)-fineK(1)
%     K4 = linspace(12, pdfK(1), 500)
%     GEV_left = gevpdf(-K4,phi,sigma,Mean)
%     nyK = [K4, nyK(2:end)]
%     nyRND = [GEV_left, nyRND(2:end)]
%     
% end
% avrg = trapz( nyK, (nyK.*nyRND) )
% snitt = [snitt avrg]
% 
% std = (nyK -avrg).^2
% std2 = trapz( nyK, (std.*nyRND) )
% avvik = [avvik std2]
% 
% M2 = std2
% 
% M3 = (nyK -avrg).^3
% M3 = trapz(nyK, (M3.*nyRND) )
% skjevhet = M3/(M2^(3/2))
% skew = [skew skjevhet]
% 
% M4 = (nyK -avrg).^4
% M4 = trapz(nyK, (M4.*nyRND) )
% kurts = M4/(M2^(2))
% kurtos = [kurtos kurts]
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
% 
% end
% 
% % xticklabels(dater)
% tel = find(~isnan(skew))
% 
% 
% %%%%%%%%%%%%%%%%%%%%
% %LOAD FUT.CSV FILE
% fid = fopen('fut.csv')
% fut = textscan(fid,'%s %s %f %f %f %f %f %f', 'Delimiter',',')
% fclose(fid)
% 
% 
% %PREPARE DATES FORMAT
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
%     snittE', avvikE', skewE', kurtosE']
% 
% nmom = sortrows(moments,1)
% 
% plot(nmom(:,1), (nmom(:,3)./T0).^0.5,nmom(:,1), nmom(:,7))
% plot(nmom(:,1), nmom(:,2),nmom(:,1), nmom(:,6))
% 
% plot(nmom(:,1), nmom(:,4),nmom(:,1), nmom(:,8))
% plot(nmom(:,1), nmom(:,5),nmom(:,1), nmom(:,9))
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