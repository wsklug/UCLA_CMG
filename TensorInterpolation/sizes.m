% This function returns the first N sizes of A, if N is a vector.  If N is a
% scalar then it returns the size of the Nth dimension.  Sometimes it is
% useful to know only the certain sizes of an array.
%
% SYNTAX: D=sizes(A,N)
%
% Example:
%
% A=rand(2,4,6,8);
%
% D=sizes(A,1:2)
% >> 2 4
%
% D=sizes(A,[1 3])
% >> 2 6
%
% DBE 2005-03-07

function D=sizes(A,N)

if length(N)==1
  D=size(A,N);
elseif isempty(N)
  D=[];
else
  for k=1:length(N)
   D(k)=size(A,N(k));
  end
%   D=D(N);
end  

return