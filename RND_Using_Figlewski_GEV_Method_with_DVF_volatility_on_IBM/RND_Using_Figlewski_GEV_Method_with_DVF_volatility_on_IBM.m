% Fits Generalized Extreme-Value tails to the RND, as per the method of
% Stephen Figlewski "Estimating the Implied Risk Neutral Density for the
% U.S. Market Portfolio".
% Fits a quadratic deterministic volatility function (DVF) to the IVs
% IV(K) = A + B*K + C*K^2;
% Fits market calls to the fitted Implied vol IV(K).
% Finds the RND from the second-order derivative of call prices with
% respect to strikes.
% Fills in the tails using the GEV distribution.
clc; clear;

% Read in IBM calls and puts to estimate the risk free and dividend rates.
Call = [...
  100.0000   48.3300   31.9250;
  105.0000   43.6500   27.1250;
  110.0000   39.4200   22.4250;
  115.0000   35.0200   17.8000;
  120.0000   31.4200   13.4500;
  125.0000   28.0200    9.4250;
  130.0000   25.1800    5.9750;
  135.0000   23.0700    3.3500;
  140.0000   21.7900    1.6700;
  145.0000   20.9700    0.7350];
Put = [...
  100.0000   42.6700    0.5850;
  105.0000   39.1000    0.7850;
  110.0000   35.3400    1.0400;
  115.0000   31.9300    1.4400;
  120.0000   28.7100    2.0550;
  125.0000   25.8600    3.0500;
  130.0000   23.2400    4.6000;
  135.0000   21.0600    6.9750;
  140.0000   19.3100   10.3000;
  145.0000   16.9500   14.3750];
K = Put(:,1);
P = Put(:,3);
C = Call(:,3);

% Spot price and Maturity.
S0 = 131.27;
T = 64/365.25;

% Risk Free rate and Dividend yield using David Shimko's approach.
[r q] = Find_RiskFree_and_Dividends_Shimko(C,P,K,S0,T);
clear K Put Call P C

% Read in call data to construct the RND.
data = [...
  100.0000   48.3300   31.9250;
  105.0000   43.6500   27.1250;
  110.0000   39.4200   22.4250;
  115.0000   35.0200   17.8000;
  120.0000   31.4200   13.4500;
  125.0000   28.0200    9.4250;
  130.0000   25.1800    5.9750;
  135.0000   23.0700    3.3500;
  140.0000   21.7900    1.6700;
  145.0000   20.9700    0.7350;
  150.0000   20.7700    0.3100];

% Extract the strikes (K) implied vols (IV) market call prices (MK).
K  = data(:,1);
IV = data(:,2)./100;
MK = data(:,3);

% Validate the implied vols.
for i=1:length(K)
 	IVcheck(i) = fminsearch(@(v) BSIV(v,MK(i),S0,K(i),r,q,T,'C'),.65);
end

[IV IVcheck']
clear IVcheck

%% Find the DVF parameters.
beta = regstats(IV, K, 'purequadratic');
yhat = beta.yhat;
B = beta.beta;
clear yhat

% Create a continuum of strikes.
dK = .1;
K2 = [min(K):dK:max(K)+2];

for i=1:length(K2);
	fVol(i) = B(1) + B(2)*K2(i) + B(3)*K2(i)^2;
end
clear B

% Plot the results.
plot(K2,fVol,K,IV,'o')
legend('DVF Implied Vols', 'Market Implied Vol')
xlabel('Strike Price')
ylabel('Implied Volatility')
pause

% Price the calls using the IV fitted from the quadratic DVF.
for i=1:length(K2);
    Call(i) = BSCall(fVol(i),S0,K2(i),r,q,T,'C');
end

% Plot the results.
plot(K2, Call, K, MK, 'o')
legend('Black-Scholes-DVF call prices', 'Market call prices')
xlabel('Strike Price')
ylabel('Call Price')
pause

% Take first derivative of market calls wrt strike.
for i=1:length(K2)-1;
    dk = K2(i+1) - K2(i);
    dc = Call(i+1) - Call(i);
    dcdk(i) = dc/dk;
end

% Take second derivative of market calls wrt strike to get the RND.
for i=1:length(dcdk)-1;
    dk = K2(i+1) - K2(i);
    dc = dcdk(i+1) - dcdk(i);
    RND(i) = dc/dk;
end

plot(K2(1:end-2), RND)
legend('DVF Risk Neutral Density')
xlabel('Stock price at maturity, S(T)');
ylabel('Risk Neutral Density')
pause

%% Integrate to check whether area = 1 (it won't).
AreaIntegralRND = trapz(RND)*dk

%% Check whether the RND recovers the call prices (it won't!).
K22 = K2(3:end);

for i=1:length(K);
    Kstar = K(i);
    I = find(K22<Kstar);
    K22(I) = 0;
    integrand = max(K22 - Kstar,0).*RND;
    RecoveredCallDVF(i) = trapz(integrand).*dk;
end;
RecoveredCallPrices = [RecoveredCallDVF' MK]

%% Fit the right tail with GEV distribution.
a1R = 0.995;
a0R = 0.95;
Ka1R = quantile(K2, a1R);
Ka0R = quantile(K2, a0R);
fKa1R = interp1(K2(1:length(RND)), RND, Ka1R);
fKa0R = interp1(K2(1:length(RND)), RND, Ka0R);

% Find the GEV parameters for the right tail.
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
Error = [checkCDF_a0 checkPDF_a0 checkPDF_a1]

% Find the GEV pdf values along the right tail.
K3 = [min(K)-10:dK:max(K)+50];
for i=1:length(K3);
	GEV_right(i) = PDF_GEV(K3(i),Mean,sigma,phi);
	if isreal(GEV_right(i)) & GEV_right(i)>0
		I(i) = i;
	else
		I(i) = 0;
	end
end
I = I(find(I~=0));

% Plot the right tail.
plot(K2(1:end-2), RND, K3(I), GEV_right(I), ':', Ka1R, fKa1R,'o',Ka0R, fKa0R,'o')
legend('DVF Risk Neutral Density', 'GEV for right tail', 'a1R anchor', 'a0R anchor')
xlabel('Stock price at maturity, S(T)');
ylabel('Risk Neutral Density')
pause

% Fit the left tail with the GEV distribution
a1L = 0.01;
a0L = 0.03;

Ka1L = quantile(K2, a1L);
Ka0L = quantile(K2, a0L);
fKa1L = interp1(K2(1:length(RND)), RND, Ka1L);
fKa0L = interp1(K2(1:length(RND)), RND, Ka0L);

% Find the GEV parameters for the left tail
% start = betaR;
options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e6, 'TolX', 1e-10, 'TolFun', 1e-10);
[betaL feval] = fminsearch(@(b) FitGEVLeftTail(b,a0L,a1L,Ka0L,Ka1L,fKa1L,fKa0L), start, options);
Mean  = betaL(1);
sigma = betaL(2);
phi   = betaL(3);

% Check how close the optimization conditions are
checkCDF_a0 = CDF_GEV(-Ka0L,Mean,sigma,phi) - (1-a0L);
checkPDF_a0 = PDF_GEV(-Ka0L,Mean,sigma,phi) - fKa0L;
checkPDF_a1 = PDF_GEV(-Ka1L,Mean,sigma,phi) - fKa1L;
Error = [checkCDF_a0 checkPDF_a0 checkPDF_a1]

% Find the GEV pdf values along the left tail.
K4 = [min(K)-50:dK:max(K)];
for i=1:length(K4);
	GEV_left(i) = PDF_GEV(-K4(i),Mean,sigma,phi);
end

% Plot the left tail.
plot(K2(1:end-2), RND, K4, GEV_left, ':', Ka1L, fKa1L,'o', Ka0L, fKa0L,'o')
legend('DVF Risk Neutral Density', 'GEV for left tail', 'a1 anchor', 'a0 anchor')
xlabel('Stock price at maturity, S(T)');
ylabel('Risk Neutral Density')
pause

% Plot the RND and both GEV pdfs
plot(K2(1:end-2), RND, ...
	 K3(I), GEV_right(I), ':',...
	 K4, GEV_left, ':',...
	 Ka1R, fKa1R,'o',...
	 Ka0R, fKa0R,'o',...
	 Ka1L, fKa1L,'o', ...
	 Ka0L, fKa0L,'o')
legend('DVF Risk Neutral Density', ...
	   'GEV for right tail', ...
	   'GEV left tail',...
       'a1R anchor',...
	   'a0R anchor',...
       'a1L anchor',...
	   'a0L anchor')
pause

%% Select only needed parts of the tails
Body  = [K2(1:end-2)' RND'];
Left  = [K4', GEV_left'];
Right = [K3(I)', GEV_right(I)'];

% Select left tail.  All observations less than min(K)
LeftTail = Left(find(Left(:,1)<min(K2)),:);

% Select right tail.  All observations greater than max(K)
RightTail = Right(find(Right(:,1)>max(K2)),:);

% Final RND.  Left and right tails grafted to body
Main = [LeftTail; Body; RightTail];

% Separate out Strikes and RND
Strike  = Main(:,1);
Density = Main(:,2);

plot(Strike, Density)
legend('Final Risk Neutral Density')
xlabel('Stock price at maturity, S(T)');
ylabel('Risk Neutral Density')
pause

% Integrate the RND
UpdatedAreaRND = trapz(Density).*dK

% Force the area under the RND to be one and check
Density = Density./UpdatedAreaRND;
Check = trapz(Density).*dK

%% Check whether the updated RND recovers the call prices
for i=1:length(K);
    Kstar = K(i);
    I = find(Strike<Kstar);
	NewStrike = Strike;
	NewStrike(I) = 0;
    integrand = max(NewStrike - Kstar,0).*Density;
    RecoveredCallDVFUpdated(i) = exp(-r*T)*trapz(integrand).*dk;
end;

RecoveredCallPrices = [MK RecoveredCallDVF' RecoveredCallDVFUpdated']

for i=1:length(Strike);
    Kstar(i) = Strike(i);
    I = find(Strike<Kstar(i));
	NewStrike = Strike;
	NewStrike(I) = 0;
    integrand = max(NewStrike - Kstar(i),0).*Density;
    NewCall(i) = exp(-r*T)*trapz(integrand).*dk;
end

% Calculate Implied Vol only for Call prices > $0.05
J = find(NewCall>0.05);
NewCall = NewCall(J);
Kstar = Kstar(J);
Strike = Strike(J);

for i=1:length(J)
	NewIV(i) = fminsearch(@(v) BSIV(v,NewCall(i),S0,Kstar(i),r,q,T,'C'),.35);
end;

plot(Strike, NewCall, K, MK,'rx-')
legend('Call Prices derived from RND', 'Original Call prices')
xlabel('Strike Price');
ylabel('Call Price')
pause

plot(Strike, NewIV, K, IV, 'rx')
legend('IV derived from RND', 'Original IV')
xlabel('Strike Price');
ylabel('Implied Volatility')

