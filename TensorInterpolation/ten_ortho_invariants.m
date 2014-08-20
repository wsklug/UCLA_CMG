function [K1,K2,K3,R1,R2,R3]=ten_ortho_invariants(T,E2,E3)
% This function calculates the six invariants that are members of the two
% orthogonal tensor invariant sets. Ki are the cylindrical invariants that
% include trace, norm_dev, and mode.  Ri are the spherical invariantes that
% include norm, FA, and mode.
%
% SYNTAX:  [K1,K2,K3,R1,R2,R3]=ten_ortho_invariants(T);
%
% INPUTS:  T - Tensor array with dimensions [3x3xLxMxN]
%
% OUTPUTS: K1 - Tensor trace array with dimensionality [LxMxN]
%          K2 - Tensor norm of the tensor deviatoric array with dimensionality [LxMxN]
%          K3 - Tensor mode array with dimensionality [LxMxN]
%          R1 - Tensor norm array with dimensionality [LxMxN]
%          R2 - Tensor FA array with dimensionality [LxMxN]
%          R3 - Tensor mode array with dimensionality [LxMxN]
%
% DBE 2005.04.25
% DBE 2008.03.14 Modified to support eigen value inputs...
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

if nargin==3
  tmp=T; clear T;
  T=zeros(3,3,size(tmp,1),size(tmp,2),size(tmp,3));
  T(1,1,:,:,:)=tmp;
  T(2,2,:,:,:)=E2;
  T(3,3,:,:,:)=E3;  
end

cout=0;  % Print messages...

devA=ten_dev(T);   % Deviatoric Tensor

K1=ten_trace(T);              % Trace
  if cout, fprintf('Done with TRACE(A).\n'); end
K2=ten_norm(devA);             % Norm(deviatoric)
  if cout, fprintf('Done with NORM(DEV(A)).\n'); end
K3=ten_mode(devA,K2);         % Mode
  if cout, fprintf('Done with MODE(A).\n'); end

R1=ten_norm(T);                % Norm
  if cout, fprintf('Done with NORM(A).\n'); end
R2=sqrt(3/2)*K2./R1;                    % FA
  if cout, fprintf('Done with FA(A).\n'); end
R3=K3;                        % Mode

% % Correct "bad" values
% if any(K1(:)<0)
%   warning('Found K1<0...correcting to equal zero.');
%   K1(K1<0)=0;
% end
% 
% if any(isnan(K2))
%   warning('Found NaN K2 values...correcting equal to zero.');
%   K2(isnan(K2))=0;
%   R2=K2./R1;
% end
% 
% if any(isnan(K3))
%   warning('Found NaN K3 and R3 values...correcting equal to zero.');
%   K3(isnan(K3))=0;
%   R3=K3;
% end
% 
% if any(isnan(R2))
%   warning('Found NaN R2 (FA) values...correcting equal to zero.');
%   R2(isnan(R2))=0;
% end

return

% ten_norm(ten_dev(T))

% DD=zeros(sizes(D,[1 2 4 3 5]));
% MM=zeros(sizes(M,[2 1 3]));
% trA=zeros(sizes(D,[4 3 5]));
% normA=zeros(sizes(D,[4 3 5]));
% normdevA=zeros(sizes(D,[4 3 5]));
% modeA=zeros(sizes(D,[4 3 5]));
% 
% for slc=1:size(D,5)
%   for j=1:3
%     for k=1:3
%       DD(j,k,:,:,slc)=rot90(squeeze(D(j,k,:,:,slc)),-1);
%     end
%   end
% end

% % Calculate the eigensystem
% [Evct,Eval]=eig_sorted(DD);

% 
% % Generate invariant maps...
% FA=ten_norm(ten_dev(DD))./ten_norm(DD);
% for slc=1:size(D,5)
%   for j=1:size(DD,3)
%     for k=1:size(DD,4)
%       trA(j,k,slc)=trace(DD(:,:,j,k,slc));
%       normA(j,k,slc)=norm(squeeze(DD(:,:,j,k,slc)),'fro');
%       normdevA(j,k,slc)=norm(ten_dev(squeeze(DD(:,:,j,k,slc))),'fro');
%       modeA(j,k,slc)=ten_mode(squeeze(DD(:,:,j,k,slc)));
%     end
%   end
%   MM(:,:,slc)=imfill(rot90(M(:,:,slc),-1),'holes');
% end
% M=MM; clear MM;
% 
% % Correct some "bad" values in the maps...
% trA(trA<0)=0;
% FA(isnan(FA))=0;
% modeA(isnan(modeA))=0;
% normdevA(isnan(normdevA))=0;
