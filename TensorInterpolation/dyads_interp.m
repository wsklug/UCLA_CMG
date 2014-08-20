function R = dyads_interp(EV,W,ED)%
% DT-MRI tensor orientation interpolation using dyadic tensor
%
% Syntax: R = dyads_interp(EV,W,ED)
%
% Inputs:
%   EV - eigenvectors for 1D-, bi- or tri-linear (3x3x 2,4 or 8)
%   W - weights (1,2 or 3)
%   ED - interpolated eigenvalues (3x3: diagonal matrix)
%
% Output:
%   R - interpolated eigenvector (3x3)
%
% Please cite the following when using this function:
%   Gahm JK and Ennis DB. Dyadic tensor-based interpolation of tensor orientation: application to car- diac DT-MRI. STACOM, MICCAI Workshop 2013, LNCS 8330, pp. 135?142, 2014.
%
% Written by JK Gahm, UCLA. 01/17/2013.
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

n = length(W);
M = zeros(3,3);

for i=1:3
  ev = squeeze(EV(:,i,:));
  if n==1
    D = (1-W)*ev(:,1)*ev(:,1)'+W*ev(:,2)*ev(:,2)';
  elseif n==2
    D = (1-W(1))*(1-W(2))*ev(:,1)*ev(:,1)'+...
        W(1)*(1-W(2))*ev(:,2)*ev(:,2)'+...
        (1-W(1))*W(2)*ev(:,3)*ev(:,3)'+...
        W(1)*W(2)*ev(:,4)*ev(:,4)';
  elseif n==3
    D = (1-W(1))*(1-W(2))*(1-W(3))*ev(:,1)*ev(:,1)'+...
        W(1)*(1-W(2))*(1-W(3))*ev(:,2)*ev(:,2)'+...
        (1-W(1))*W(2)*(1-W(3))*ev(:,3)*ev(:,3)'+...
        W(1)*W(2)*(1-W(3))*ev(:,4)*ev(:,4)'+...
        (1-W(1))*(1-W(2))*W(3)*ev(:,5)*ev(:,5)'+...
        W(1)*(1-W(2))*W(3)*ev(:,6)*ev(:,6)'+...
        (1-W(1))*W(2)*W(3)*ev(:,7)*ev(:,7)'+...
        W(1)*W(2)*W(3)*ev(:,8)*ev(:,8)';
  end
  
  % find the largest eigenvector
  [v,d] = teneig(D);
  M(:,i) = v(:,1);
end

% closest orthogonal matrix
if nargin<3
  ED = eye(3);
end
N = M*ED;
[U,S,V] = svd(N);
R = U*V';

if det(R)<0
  R = -R;
end



