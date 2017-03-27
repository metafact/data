
Opt = cell(1,3)
Fut = cell(1,8)
dinfo = dir('*.CSV');
for K = 1 : length(dinfo)
  thisfilename = dinfo(K).name; %just the name
  if length(dinfo(K).name) > 12 
      fileID = fopen(thisfilename)
      O = textscan(fileID,'%s %s %f', 'Delimiter',',')
      fclose(fileID)
        for i = 1:3
            Opt{i} = [Opt{i}; O{i}]
        end 
  else
      fileID = fopen(thisfilename)
      F = textscan(fileID,'%s %s %f %f %f %f %f %f', 'Delimiter',',')
      fclose(fileID)
        for i = 1:8
            Fut{i} = [Fut{i}; F{i}]
        end 
  end
end

tic
parpool('local',4) %%%
dinfo = dir('*.CSV');
Opt = cell(1,3, length(dinfo))
O = cell(1,3)
parfor K = 1 : length(dinfo)
  thisfilename = dinfo(K).name; %just the name
  if length(dinfo(K).name) > 12 
      fileID = fopen(thisfilename)
      Opt(:,:,K) = textscan(fileID,'%s %s %f', 'Delimiter',',')
      fclose(fileID)
            
  end
end
delete(gcp('nocreate'))
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
parpool('local',4) %%%
dinfo = dir('*.CSV');
Opt = cell(1,3, length(dinfo))
O = cell(1,3)
parfor K = 1 : length(dinfo)
  thisfilename = dinfo(K).name; %just the name
  if length(dinfo(K).name) > 12 
      fileID = fopen(thisfilename)
      Opt(:,:,K) = textscan(fileID,'%s %s %f', 'Delimiter',',')
      fclose(fileID)
            
  end
end
delete(gcp('nocreate'))
toc


%21.03.2017
E = cell(1,3)
dinfo = dir('*.CSV');
for K = 1 : length(dinfo)
  thisfilename = dinfo(K).name;  %just the name
  fileID = fopen(thisfilename)
  C = textscan(fileID,'%s %s %f', 'Delimiter',',')
  fclose(fileID)
    for i = 1:3
        E{i} = [E{i}; C{i}]
    end 
end

idx = find(ismember(E{1}, 'CL2009Z5100C'))
E{1}(idx)
op = strmatch('CL2009Z51', E{1})
 op = strmatch(F{2}(800), E{2})
  pp = strmatch(F{1}(1), E{1})
C = intersect(op,pp)

% FOR CELLS, USEFUL FUNCTIONS
-cellfun
-accumarray

optC = strcat(opt{1}, opt{2})


for K = 1 : length(dinfo)
thisfilename = dinfo(K).name;
dinfo = dir('*.CSV');
%%%%%%%%%%%%!!!!!!!!!START
 fid = fopen('new.csv')
      opt = textscan(fid,'%s %s %f', 'Delimiter',',')
 fclose(fid)
 
 fid = fopen('fut.csv')
      fut = textscan(fid,'%s %s %f %f %f %f %f %f', 'Delimiter',',')
 fclose(fid)
opReg = regexprep(opt{1}(1:end-1), '.{5}$','')
optC = strcat(opReg, opt{2}); % Because there is some sign at the 'end' index
futC = strcat(fut{1}(1:end-1), fut{2});


[~, loc] = ismember(optC,futC) %R is concatened cell array of options and W is for futures, this function assignes index from W array 
%for each element from array R
lo = loc(loc ~= 0); %deletes elements with zeros
%lo = find(loc ~= 0) % the same as previuos command. It obtains indicies where 'loc' elements are not zero

%Next three lines are for preparing R array
z = find(loc == 0)
i = 1:length(optC);
i(z) = []; %deletes indeces where elements are zeros


%TO CHECK IF BOTH ARE EQUAL
optCor = optC(i);
futCor = futC(lo);
isequal(optCor,futCor) %1 if they are equal

%Now get Strikes
A = cell2mat(opt{1}(i));
K = A(:,8:11);
K = str2num(K)/100;

%Ma = cell2mat((opt{1}(1:end-1)));
%Ma = Ma(:,3:7);
Ma = opt{1}(i)
Mat = regexprep(Ma, 'Z', '/11/19');
Mat1 = regexprep(Mat, 'K', '/04/19');
Mat2 = regexprep(Mat1, 'G', '/01/19');
Mat3 = regexprep(Mat2, 'H', '/02/19');
Mat4 = regexprep(Mat3, 'J', '/03/19');
Mat5 = regexprep(Mat4, 'M', '/05/19');
Mat6 = regexprep(Mat5, 'N', '/06/19');
Mat7 = regexprep(Mat6, 'Q', '/07/19');
Mat8 = regexprep(Mat7, 'U', '/08/19');
Mat9 = regexprep(Mat8, 'V', '/09/19');
Mat10 = regexprep(Mat9, 'X', '/10/19');
Mat11 = regexprep(Mat10, 'F', '/12/19');

j = contains(Mat11, '/12/19')
last = cell2mat(Mat11);
last = last(:,3:12);
t  = last(j,:)
ar = str2num(t(:,1:4))
arM1 = ar-1
arM1 = num2str(arM1)
dat = strcat(arM1, t(:,5:end))
last(j,:) = dat

%MAKE A FINAL TABLE
 T = table(fut{6}(lo), K, opt{3}(i), opt{2}(i), last)
 save('table.mat', 'T') 
 writetable(T, 'table.csv')

cellfun(@(d) addtodate(d, -1, 'year'), x)
 %{
 SOME MANIPULATION WITH ONE YEAR SUBSTRACTION
 
 t  = datenum('2009/11/19') ->get date to number
t= addtodate(t, -1, 'year') -> substracts one year 
 datestr(t) -> returns to string format
 
 %}






S = F{6}(lo) -K
%%%%%%%%%

'


for i = 1:length(A)
	switch A(i)
		case 'Z'
			D (i) = 
if 

%2. Aproxaimte with barone Adesi and Whaley

%3. Read for the all CL/CB options and Futures


E{4} = [] % creates 4th cell array as a column to original E cell array
C = cell(11996,1) % creates cell array with 11996X1 dims
C(:) = {'String'}  %adds STRING elements to cell array, so ellemens are all the same

v = 0.03
m=repmat(v,10889,1) %creates 10889x1 array with the same v element repeating everywhere

[row, col] = find(isnan(YourMatrix)) % Look for NaN elements in array

plot3(x,y,z, '.', 'MarkerSize', 50) %Plots 3D graph without any functional specifications, '.' is a point plotting 
grid on								% and MarkerSize is usual size of a point

%fig2plotly() : COOL 3D PLOT FOR MATLAB


%%% BASH COMMAND FOR REPLACING SOME TEXT
gc myFile.txt) -replace 'foo', 'bar' | Out-File myFile.txt"
%%%

%%%%%%%%%%%%{
dinfo = dir('*.CSV');
tic
for K = 1 : length(dinfo)
    thisfilename = dinfo(K).name;
     fid = fopen(thisfilename)
      opt = textscan(fid,'%s %s %f', 'Delimiter',',')
     fclose(fid)

     opReg = regexprep(opt{1}, '.{5}$','')
     optC = strcat(opReg, opt{2}); % Because there is some sign at the 'end' index
     futC = strcat(fut{1}(1:end-1), fut{2});
     [~, loc] = ismember(optC,futC)
     lo = loc(loc ~= 0);

     z = find(loc == 0)
     i = 1:length(optC);
     i(z) = [];

     optCor = optC(i);
    futCor = futC(lo);
    isequal(optCor,futCor)
    
    A = cell2mat(opt{1}(i));
    K = A(:,8:11);
    K = str2num(K)/100;

    %Ma = cell2mat((opt{1}(1:end-1)));
    %Ma = Ma(:,3:7);
   
 t = cell2mat(opt{1}(1));
    n = t(7)
    
    switch n
        case 'Z'
        tt = strcat(t(3:6), '/11/19')
        case 'K'
        tt = strcat(t(3:6), '/04/19')
        case 'G'
        tt = strcat(t(3:6), '/01/19')
        case 'H'
        tt = strcat(t(3:6), '/02/19')
        case 'J'
        tt = strcat(t(3:6),  '/03/19')
        case 'M'
        tt = strcat(t(3:6), '/05/19')
        case 'N'
        tt = strcat(t(3:6), '/06/19')
        case 'Q'
        tt = strcat(t(3:6), '/07/19')
         case 'U'
        tt = strcat(t(3:6), '/08/19')
         case 'V'
        tt = strcat(t(3:6), '/09/19')
        case 'X'
        tt = strcat(t(3:6), '/10/19')
        case 'F'
        tt = strcat(num2str(str2num(t(3:6))-1), '/12/19')
    end
   rate = repmat(0.04, size(i,2),1)
   T = table(fut{6}(lo), K, opt{3}(i), cell2mat(opt{2}(i)), repmat(tt, size(i,2),1),rate)
    vol = {};
   for j = 1:size(i,2)
       RateSpec = intenvset('ValuationDate', T.Var4(j,:), 'StartDates', T.Var4(j,:),...
'EndDates', T.Var5(j,:), 'Rates', T.rate(j), 'Compounding', -1, 'Basis', 1)
        StockSpec = stockspec(NaN, T.Var1(j), {'continuous'}, 0.1)
        OptSpec = {'call'};
        ImpVol =  impvbybaw(RateSpec, StockSpec, T.Var4(j,:), T.Var5(j,:), OptSpec,...
T.K(j), T.Var3(j))
        vol = [vol, ImpVol]
   end 
   Volat = cell2mat(vol)
   Volat
 toc   
end
    
%%%%%%%%%%%%%%%}
