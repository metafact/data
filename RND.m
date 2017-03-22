
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


[~, loc] = ismember(R,W) %R is concatened cell array of options and W is for futures, this function assignes index from W array 
%for each element from array R
lo = loc(loc ~= 0) %deletes elements with zeros
lo = find(loc ~= 0) % the same as previuos command. It obtains indicies where 'loc' elements are not zero

%Next three lines are for preparing R array
z = find(loc == 0)
i = 1:length(R)
i(z) = [] %deletes indeces where elements are zeros



Rr = R(i)
Ww = W(lo)

isequal(Rr,Ww) = 1 %1 if they are equal

%Now get Strikes
A = cell2mat(E{1}(i))
K = A(:,8:11)
K = str2num(K)/100

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