function y = PDF_GEV(S,Mean,sigma,phi)

z = (S-Mean)/sigma;

y = CDF_GEV(S,Mean,sigma,phi) * (1+phi*z)^(-1/phi-1) / sigma;

