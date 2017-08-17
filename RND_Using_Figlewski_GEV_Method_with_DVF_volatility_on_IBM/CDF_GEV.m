function y = CDF_GEV(S,Mean,sigma,phi)

z = (S-Mean)/sigma;

y = exp(-(1+phi*z)^(-1/phi));

