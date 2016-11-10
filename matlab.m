d = fetch(Connect, {'SP5SOGE','SP5SIOG'}, 'PI', '12/31/1980', '12/31/2015', 'Q' ) %SP5SOGE -Exploration and Production, SP5IOG -Integrated

class(a) %-class of an object a
d(1).PI  %-Closing prices for first index

cd M: %change to disk M

help <command> % to get descrition of the command

xlswrite('test.xls', N, 'Quarterly') %writes to excel file, but first Create file to write into, N-matrix or table to write, 
%'Quarterly' writes to this tab

fid = fopen('test.csv', 'w')
fprintf(fid, '%s,', data.DATE)
