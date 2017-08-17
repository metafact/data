% Check FMINSEARC -> p0 and p1

snitt = []
avvik = []
skew = []
kurtos = []

datoer = []
ttm = []

snittTS = []
avvikTS = []
skewTS = []
kurtosTS = []

datoerTS = []
ttmTS = []

% teller = 0
neg = find(Uniq.Vols < 0)
Uniq.MatCell = cellstr(Uniq.Maturity); %Adding extra column which is converted Maturity dates
% m = find(Uniq.MatRatio == 30/365); %Selecting only one month Maturity
% Unik = sortrows(unique(Uniq.MatCell(m))) %Now sorting the unique MATURITY rows ACSENDING
Unik = sortrows(unique(Uniq.MatCell))
for i = 1:length(Unik) %For each unique 30 days Maturity date
    j = find(contains(Uniq.MatCell, Unik(i)) & Uniq.MatRatio <= 30/365) %Indeces of DATES for one specific 30 MATURITY
    nt = Uniq(j,:); %Table with ontly <= MatRatio
    
    nt.DateNum = datenum(nt.Date); %Convert to date numerics for sorting, b/c sorting for CELL format is not working
    nt = sortrows(nt, 10)%New table with each specific MATURITY date
    
    datUniq = unique(nt.DateNum) %get unique dates in NUMERIC format
    formatOut = 23 %format for mm/dd/yyyy
    datUniq = datestr(sortrows(datUniq), formatOut) %Unique dates for new table with one one same MATURITY
    
    
    %parfor
    parfor m = 1:size(datUniq,1)
        mj = find(contains(nt.Date, datUniq(m,:))) %INDECES for only one specific DATE and one specific MATURITY or <=30
        if size(mj) < 5
            %             continue
            snitt(m) = NaN
            avvik(m) = NaN
            skew(m) = NaN
            kurtos(m) = NaN
            
            datoer(m) = NaN
            ttm(m) = NaN
        else
            [avrg, std2, skjevhet, kurts] = RiNeDe(nt(mj,:))
            
            snitt(m) = avrg
            avvik(m) = std2
            skew(m) = skjevhet
            kurtos(m) = kurts
            
            datoer(m) = datenum(datUniq(m,:))
            ttm(m) = (datoer(m) - datenum( Unik(i) )) / 365
        end
    end
    snittTS = [snittTS snitt]
    avvikTS = [avvikTS avvik]
    skewTS = [skewTS skew]
    kurtosTS = [kurtosTS kurtos]
    
    datoerTS = [datoerTS datoer]
    ttmTS = [ttmTS ttm]
    
    snitt = []
    avvik = []
    skew = []
    kurtos = []
    
    datoer = []
    ttm = []
end

tab = table(datoerTS', snittTS', avvikTS', skewTS', kurtosTS', ttmTS')
na = find(isnan(tab.Var2))
tab(na,:) = []
srt = sortrows(tab,1)
xtickangle(90)
datetick('x',20)
srt.Properties.VariableNames{'Var1'} = 'Dates'
srt.Properties.VariableNames{'Var2'} = 'Mean'
srt.Properties.VariableNames{'Var3'} = 'Vol'
srt.Properties.VariableNames{'Var4'} = 'Skew'
srt.Properties.VariableNames{'Var5'} = 'Kurt'
srt.Properties.VariableNames{'Var6'} = 'TTM'
srt.DateStr = datestr(srt.Dates)

plot(srt.Dates, (srt.Vol.^0.5).*sqrt(252)) %Plots volatility

srt.Properties.VariableNames
srt.VolsAnnual = (srt.Vol.^0.5).*sqrt(252);
srt.DateStr = datestr(srt.Dates);
writetable(srt, 'Moments.csv')


