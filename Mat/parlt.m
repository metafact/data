
Tab.Vol = zeros(length(Tab.K),1)

tic
parfor i = 1:size(Tab,1)
    vol = impliedV(Tab.Var1(i), Tab.K(i), Tab.time(i), 0.04, Tab.Var3(i))
    Tab.Vol(i) = vol(1)
    %pris =[pris pr]
end
toc
%{
x0 = [1,1]'
options = optimoptions('lsqnonlin', 'display', 'none');
x=lsqnonlin(@objFunctions,x0,[0.01,1],[],options);


function [ rssOutput ] = objFunctions( params) 
F = 50
r = 0.03
time = 0.5
b = 0
X = 40
AmPr = 1.8

sigma = params(1)
S_star = params(2)

sigma_sqr = sigma*sigma;
time_sqrt = sqrt(time);
nn = 2*b/sigma_sqr; 
m = 2*r/sigma_sqr;  
K = 1-exp(-r*time); 
q2 = (-(nn-1)+sqrt((nn-1)^2+(4*m/K)))*0.5

e1 = (log((F/S_star)/X)+(b+0.5*sigma_sqr)*time)/(sigma*time_sqrt);
B2 =  ((1-exp((b-r)*time)*normcdf(e1))* ( (F/S_star) /q2))^q2; 
delta = AmPr -blsprice(F, X,r,time, sigma)-B2

d1 = (log(S_star/X)+(b+0.5*sigma_sqr)*time)/(sigma*time_sqrt);
    A2 =  (1-exp((b-r)*time)*normcdf(d1))* (S_star/q2); 
    fu = S_star - X -blsprice(S_star, X,r,time, sigma) - A2

    
    rssOutput = [delta; fu]; 
  end 
%}
