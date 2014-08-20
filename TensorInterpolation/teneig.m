function [V,D] = teneig(T)
% Find eigenvalues and eigenvectors of tensors
%
% Input
%  - T is a diffusion tensor (3x3xLxMxN or 3x3x(L*M*N)).
% Output
%  - V is eigenvectors (3x3xLxMxN).
%  - D is eigenvalues (3x3xLxMxN).
%
% Written by Jin Kyu Gahm, UCLA, July 2009
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

[x1,x2,L,M,N] = size(T);
numVoxels = L*M*N;

V = zeros(3,3,L,M,N);
D = zeros(3,3,L,M,N);

for i=1:numVoxels
    A = T(:,:,i);
    [VA,DA] = eig(A);
    [dA,j] = sort(diag(DA),'descend');
    D(:,:,i) = diag(dA);
    V(:,:,i) = VA(:,j);
    
    if det(V(:,:,i))<0
      V(:,:,i) = V(:,:,i)*diag([1 1 -1]);
    end
    
%     if i>1 && mod(i-1,L*M)==0
%         fprintf(['slice ',int2str(i/(L*M)),' done!\n']);
%     end
end


%% not efficient for large memory
% T = reshape(T,[x1 x2 L*M N]);
% V = zeros(3,3,L,M,N);
% D = zeros(3,3,L,M,N);
% 
% %fprintf('Computing eigenvectors\n');
% for slc=1:N
%     for i=1:L*M
%         A = T(:,:,i,slc);
%         [VA,DA] = eig(A);
%         [dA,j] = sort(diag(DA),'descend');
%         D(:,:,i,slc) = diag(dA);
%         V(:,:,i,slc) = VA(:,j);
%     end
%     %fprintf('.');
% end
% %fprintf('\n');

