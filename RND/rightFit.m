function [indeks, prob] = rightFit(fineK, newC, rf, T0, alpha, t)
    %CALCULATE INDEX FOR GIVEN ALPHA THREASHOLD FIRST
    try
        
		i = 2
		p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
		p = p * exp(rf*T0) +1
		while p < alpha  % Mojet ot desyati procentov nachatsya DAJE??
			i = i + 1
			p = (newC(i+1)-newC(i-1))/ (fineK(i+1)-fineK(i-1))
			p = p * exp(rf*T0) +1
        end
    end %try
    
    %TRY DIFFERENT CASES 
	if i == 500 & t == 1
        if alpha ~= 0.85
            alpha = 0.85
            [indeks, prob] = rightFit(fineK, newC, rf, T0, alpha, t)
           
        else
            indeks = NaN
            prob = NaN
        end
% 		alpha = 0.04
%       leftFit(fineK, newC, rf, T0, alpha)
    elseif i == 500 & t == 2
        indeks(t) = 500
        prob(t) = 0
        
    elseif i < 500  
        
        indeks(t) = i
        prob(t)= p
        t = t + 1
        %CALL RECURSION FUNCTION IF 'i' IS LOCATED QUITE HIGH IN OUR RND
        %CURVE
        if t  < 3
            [itemp, ptemp] = rightFit(fineK, newC, rf, T0, alpha + 0.03, t)
                         
        else
                return
        end
    indeks(t) = itemp(t)
    prob(t) = ptemp(t)
    end %if
    
    
end

% ze = find(~cellfun(@isempty, ks))
% for i = 1:20
%     figure(i)
%     plot(ks{ ze(i) }, RNDs{ ze(i) })
% end
