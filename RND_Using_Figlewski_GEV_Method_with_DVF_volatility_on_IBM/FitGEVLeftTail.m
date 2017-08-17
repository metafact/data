function y = FitGEVLeftTail(beta, a0, a1, Ka0, Ka1, fKa1, fKa0)

Mean  = beta(1);
sigma = beta(2);
phi   = beta(3);

y = (CDF_GEV(-Ka0,Mean,sigma,phi) - (1-a0))^2 + ...
	(PDF_GEV(-Ka0,Mean,sigma,phi) - fKa0)^2 + ...
	(PDF_GEV(-Ka1,Mean,sigma,phi) - fKa1)^2 ;

if sigma<0 | isreal(y)==0
	y = 1e100;
end



