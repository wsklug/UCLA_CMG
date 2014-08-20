% This function generates the deviatoric (anisotropic) part of the tensor A.  
%
% SYNTAX:  D=ten_dev(A)
%
% INPUTS:  A - A tensor with dimensions [Nx3x3] or [LxMxNx3x3]
%
% OUTPUTS: D - Matches INPUTS dimensionality.
%
% DBE@UCLA 02.05.2002
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

function D=ten_dev(A)

N_dims_A=ndims(A);
A_sz=sizes(A,3:N_dims_A);

if N_dims_A==5 | N_dims_A==4  % Multi- or single-slice
  A_tmp=reshape(A,[3 3 prod(A_sz)]);
elseif N_dims_A==3 | N_dims_A==2
  A_tmp=A;
end

clear A;

D=zeros(3,3,size(A_tmp,3));
for k=1:size(A_tmp,3)
  D(:,:,k)=-(1/3)*trace(A_tmp(:,:,k))*eye(3);
end

D=A_tmp+D;

if N_dims_A==5 | N_dims_A==4  % Multi- or single-slice
  D=reshape(D,[3 3 A_sz]);  
end

% if N_dims_A==5
%   M=zeros(size(A));
%   for i=1:size(A,1);
%     for j=1:size(A,2);
%       for k=1:size(A,3)
%         D(i,j,k,:,:)=squeeze(A(i,j,k,:,:))-(1/3)*trace(squeeze(A(i,j,k,:,:)))*eye(size(squeeze(A(i,j,k,:,:))));
%       end
%     end
%   end
% elseif N_dims_A==3
% 	for k=1:size(A,3)
% 		D(:,:,k)=A(:,:,k)-(1/3)*trace(A(:,:,k))*eye(size(A(:,:,1)));
% 	end
% end

return
