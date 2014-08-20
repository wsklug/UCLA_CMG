% This function computes the tensor's Frobenius norm of an array that is has
% dimensions [KxKxLxMxN].
%
% SYNTAX:  M=ten_norm(A)
%
% INPUTS:  A - An array of tensor data of dimension [KxKxLxMxN]
%
% OUTPUTS: M - An array of tensor magnitudes of dimension [LxMxN]
%
% Tensor magnitude is defined as sqrt(A:B)=sqrt(trace(A*B'))
%
% DBE 2005.01.19
% DBE 2006.08.31 Cleaned up and renamed from ten_mag to ten_norm
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

function M=ten_norm(A)

if ndims(A)==5 | ndims(A)==4 % Multi- or single-slice
  A_tmp=reshape(A,[3 3 prod(sizes(A,3:ndims(A)))]);
elseif ndims(A)==3 | ndims(A)==2
  A_tmp=A;
end

M=zeros(1,size(A_tmp,3));
for k=1:size(A_tmp,3)
  M(k)=norm(A_tmp(:,:,k),'fro');
end

if ndims(A)==5 | ndims(A)==4
  M=reshape(M,sizes(A,3:ndims(A)));  
end

return

% % M=A_tmp;
% 
% % if ndims(A)==5
% %   M=zeros(sizes(A,1:3));
% %   for i=1:size(A,1);
% %     for j=1:size(A,2);
% %       for k=1:size(A,3)
% %         M(i,j,k)=norm(squeeze(A(i,j,k,:,:)),'fro');
% %       end
% %     end
% %   end
% % elseif ndims(A)==3
% %   for k=1:size(A,3)
% %     M(:,:,k)=norm(A(:,:,k),'fro');
% %   end
% % end
% 
% 
% 
% % % for k=1:size(A,3)
% % %   tmp1=squeeze(A(:,:,k));
% % %   tmp2=squeeze(A(:,:,k))';
% % %   M=sqrt(trace(tmp1*tmp2))
% % % %   M=sqrt(trace(A(:,:,k)*A(:,:,k)'));
% % % end
% % 

