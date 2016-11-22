
Connect = datastream('DS:User1','Pass1','Datastream','http://dataworks.thomson.com/Dataworks/Enterprise/1.0')
d = fetch(Connect, {'SP5SOGE','SP5SIOG', 'S&PCOMP'}, 'PI', '12/31/1980', '12/31/2015', 'Q' ) %SP5SOGE -Exploration and Production, SP5IOG -Integrated, 
S&PCOMP -index 500 itself, USTBL3M -US Treasury Bills 3 month, not a secondary

class(a) %-class of an object a
d(1).PI  %-Closing prices for first index

cd M: %change to disk M

help <command> % to get descrition of the command

xlswrite('test.xls', N, 'Quarterly') %writes to excel file, but first Create file to write into, N-matrix or table to write, 
%'Quarterly' writes to this tab

fid = fopen('test.csv', 'w')
fprintf(fid, '%s,', data.DATE)

clipboard('copy',data)% copies data to the clipboard. 
 %If data is not a character array, it is converted using mat2str.
txt = clipboard('paste') %returns the current contents of the clipboard as a character vector. 
%If clipboard cannot convert the contents, txt is empty ('').

data = clipboard('pastespecial') 5imports the clipboard contents into an array using uiimport.

