function [r q] = Find_RiskFree_and_Dividends_Shimko(C,P,K,S0,T)

% Estimates of the risk free rate and dividend yield from David Shimko.
% Requires inputs:
% C = Call prices
% P = Put prices
% K = Strikes
% S0 = Spot Price
% T = Time to maturity
% Returns:
% r = Estimate of the risk free rate
% q = Estimate of the dividend yield

% Regression to get dividend yield and risk-free rate
B = regstats(C-P,K,'linear');
int   = B.beta(1);
slope = B.beta(2);

% Risk free rate
r = -log(-slope)/T;
if r<0
	r=0;
end

% Dividend yield
D = int/S0;
q = -log(D)/T;
if q<0
	q=0;
end

% Shimko values
int = 387.9812;
slope = -0.991646;
T = 2/12;
S0 = 390.02;

