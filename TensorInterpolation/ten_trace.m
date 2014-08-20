% This function calculates the trace of tensors in an ND-arrays (3x3xXxYxZ).  
% The trace is taken of the first two dimensions and the subsequent dimensions are
% assumed to be spatial coordinates.
%
% SYNTAX:  trA=ten_trace(T);
%
% INPUTS:  T - Array of tensors with dimensionality of [3x3xXxYxZ]
%
% OUTPUTS: trA - Array of tensor trace values with dimensionality of [XxYxZ]
%
% DBE@UCLA 2005.04.25
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

function M=ten_trace(A);

% if ndims(T)==2
%   trA=0;
% else
%   trA=zeros(sizes(T,3:ndims(T)));
% end
% 
% for j=1:size(T,3)
%   for k=1:size(T,4)
%     for l=1:size(T,5)
%       trA(j,k,l)=trace(T(:,:,j,k,l));
%     end
%   end
% end

if ndims(A)==5 | ndims(A)==4 % Multi- or single-slice
  A_tmp=reshape(A,[3 3 prod(sizes(A,3:ndims(A)))]);
elseif ndims(A)==3 | ndims(A)==2
  A_tmp=A;
end

M=zeros(1,size(A_tmp,3));
for k=1:size(A_tmp,3)
  M(k)=trace(A_tmp(:,:,k));
end

if ndims(A)==5 | ndims(A)==4
  M=reshape(M,sizes(A,3:ndims(A)));  
end


return