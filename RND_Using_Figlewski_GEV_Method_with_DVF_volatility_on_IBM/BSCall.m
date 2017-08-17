function y = BSCall(v,S,K,r,q,T,PutCall);

d1 = (log(S/K) + T*(r-q+v^2/2)) / v / sqrt(T);
d2 = d1 - v*sqrt(T);

call = S*exp(-q*T)*normcdf(d1) - exp(-r*T)*K*normcdf(d2);
put = call + K*exp(-r*T) - S*exp(-q*T);

if PutCall=='C'
	y = call;
else
	y = put;
end
