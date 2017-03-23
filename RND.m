
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


%%%%%%%%%%%%
opReg = regexprep(opt{1}(1:end-1), '.{5}$','')
optC = strcat(opReg, opt{2}) % Because there is some sign at the 'end' index
futC = strcat(fut{1}(1:end-1), fut{2})


[~, loc] = ismember(optC,futC) %R is concatened cell array of options and W is for futures, this function assignes index from W array 
%for each element from array R
lo = loc(loc ~= 0) %deletes elements with zeros
lo = find(loc ~= 0) % the same as previuos command. It obtains indicies where 'loc' elements are not zero

%Next three lines are for preparing R array
z = find(loc == 0)
i = 1:length(optC);
i(z) = []; %deletes indeces where elements are zeros



optCor = optC(i);
futCor = futC(lo);

isequal(optCor,futCor) %1 if they are equal

%Now get Strikes
A = cell2mat(optC(i))
K = A(:,8:11);
K = str2num(K)/100;

Ma = cell2mat((opt{1}(1:end-1)));
Ma = Ma(:,3:7);
 Mat = regexprep(Ma, 'Z', '/11/19');

 %{
 SOME MANIPULATION WITH ONE YEAR SUBSTRACTION
 
 t  = datenum('2009/11/19') ->get date to number
t= addtodate(t, -1, 'year') -> substracts one year 
 datestr(t) -> returns to string format
 
 %}






S = F{6}(lo) -K
%%%%%%%%%

%1. How to match for Letters of Months
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

