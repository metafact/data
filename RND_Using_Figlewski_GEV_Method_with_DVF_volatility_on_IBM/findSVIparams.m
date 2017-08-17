function y = findSVIparams(params,K,IV,S0)

a = params(1);
b = params(2);
r = params(3);
m = params(4);
s = params(5);

for i=1:length(K);
    modelVol(i) = SVI(K(i),S0,a,b,r,m,s);
    e(i) = (modelVol(i) - IV(i))^2;
end

y = sum(e);


