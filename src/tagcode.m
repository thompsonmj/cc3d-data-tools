function id = tagcode
%
% TAGCODE Create an 8-character alphanumeric ID code.
%
% Description 
%   This function creates a 8-character alphanumeric code used to tag simulation parameter sets and
%   their resulting data output and analyses.

charList = ['A':'Z' '0':'9'];
randNList = randi(numel(charList),[1,8]);
id = charList(randNList);

end