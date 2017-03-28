function [impvol] = impliedV(F, X, time, r, AmPr)
b = 0
x0 = [1,1]'
options = optimoptions('lsqnonlin', 'display', 'none');
x=lsqnonlin(@objFunctions,x0,[0.01,1],[],options);

impvol = x
function [ rssOutput ] = objFunctions( params) 


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

end

%{
function impvol = impliedV(F, X, time, r, optVal)
    b = 0
   
    
    S_star = 0
    sigma = 0.01
    
    sigma_sqr = sigma*sigma;
    time_sqrt = sqrt(time);
    nn = 2*b/sigma_sqr; 
    m = 2*r/sigma_sqr;  
    K = 1-exp(-r*time); 
    q2 = (-(nn-1)+sqrt((nn-1)^2+(4*m/K)))*0.5
    
    options = optimoptions('lsqnonlin', 'display', 'none');
    [S_star,resnorm] = lsqnonlin(@pris,S_star,1,[],options);

   
    [sigma,resnorm] = lsqnonlin(@volat,sigma,0,[],options); 
    
    
    function fu = pris( S_star )
    d1 = (log(S_star/X)+(b+0.5*sigma_sqr)*time)/(sigma*time_sqrt);
    A2 =  (1-exp((b-r)*time)*normcdf(d1))* (S_star/q2); 
    fu = S_star - X -blsprice(S_star, X,r,time, sigma) - A2
    end


%{
if F < S_star
    options = optimoptions('lsqnonlin', 'display', 'none');
    [sigma,resnorm] = lsqnonlin(@volat,sigma,[],[],options); 
end
%}
    %
    function delta = volat(S_star)

       sigma_sqr = sigma*sigma;
       %time_sqrt = sqrt(time);
       nn = 2*b/sigma_sqr; 
       m = 2*r/sigma_sqr;  
       %K = 1-exp(-r*time); 
       q2 = (-(nn-1)+sqrt((nn-1)^2+(4*m/K)))*0.5

       e1 = (log((F/S_star)/X)+(b+0.5*sigma_sqr)*time)/(sigma*time_sqrt);
        B2 =  ((1-exp((b-r)*time)*normcdf(e1))* ( (F/S_star) /q2))^q2; 
        %delta = F - X
        %
        if F <S_star
            delta = optVal -blsprice(F, X,r,time, sigma)-B2
        else
            delta = optVal-F + X
        
        end
       
         %delta = optVal -blsprice(S_star, X,r,time, sigma)-B2
    end
    impvol = sigma
end
%}
