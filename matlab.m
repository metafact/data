%{

This change has been incorporated into the documentation in Release 2011a (R2011a). For previous releases, read below for any additional information:

To get the static data from a list, you can type :

x = fetch(datastream_connection,'LFRCAC40~REP','P')
or

x = fetch(datastream_connection,'LFRCAC40~REP~=P')
The REP is needed to denote a static report. This is a Datastream specific syntax.

NAME, SECD, and ISIN are some examples of the static fields and need to be requested

with the '~REP' flag:

data = fetch(datastream_connection,{'IBM~REP'}, {'NAME'});
%}

% EDIT <FUNCTION NAME> => to see the source code of the function


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


imp =xlsread('Bok1.xlsx') % Read excell file
fac = imp(:,1:21) %partione
stepwisefit(fac,imp(:,24), 'penter',0.05,'premove',0.10) %premove-factors with p value bigger than 'premove' are removed.
%penter -variables with p value lower than 'penter' are included into regression 
