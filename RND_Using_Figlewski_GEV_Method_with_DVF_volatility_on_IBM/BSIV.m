function y = BSIV(v,Price,S,K,r,q,T,PutCall);

d1 = (log(S/K) + T*(r-q+v^2/2)) / v / sqrt(T);
d2 = d1 - v*sqrt(T);

call = S*exp(-q*T)*normcdf(d1) - exp(-r*T)*K*normcdf(d2);
put = call + K*exp(-r*T) - S*exp(-q*T);

if PutCall=='C'
	BSPrice = call;
else
	BSPrice = put;
end

y = (Price - BSPrice)^2;

if v<0
	y = 1e10;
end
