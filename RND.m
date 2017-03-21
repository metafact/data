
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


%Next three lines are for preparing R array
z = loc(loc == 0)
i = 1:length(R)
i(z) = [] %deletes indeces where elements are zeros

lo = find(loc ~= 0) % the same as two previuos commands. It obtains indicies where 'loc' elements are not zero

Rr = R(i)
Ww = W(lo)

isequal(Rr,Ww) = 1 %1 if they are equal

%Now get Strikes
A = cell2mat(E{1}(i))
K = A(:,8:11)
K = str2num(K)/100

S = F{6}(lo) -K


