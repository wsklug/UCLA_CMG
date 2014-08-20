% This function computes the TENSOR determinant of an array that is has
% dimensions [3x3xN] or [3 3 X Y Z].
%
% SYNTAX: M=ten_det(A)
%
% DBE 2005/07/13

function M=ten_det(A)

if ndims(A)==5 | ndims(A)==4 % Multi- or single-slice
  A_tmp=reshape(A,[3 3 prod(sizes(A,3:ndims(A)))]);
elseif ndims(A)==3 | ndims(A)==2
  A_tmp=A;
end

M=zeros(1,size(A_tmp,3));
for k=1:size(A_tmp,3)
%   M(k)=norm(A_tmp(:,:,k),'fro');
  M(k)=det(A_tmp(:,:,k));
end

if ndims(A)==5 | ndims(A)==4
  M=reshape(M,sizes(A,3:ndims(A)));  
end

return