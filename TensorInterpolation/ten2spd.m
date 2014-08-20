function Dnew = ten2spd(D,type)
% Convert non-SPD (symmetric positive definite) tensors into SPD tensors
%
% Syntax: Dnew = ten2spd(D,type)

% Inputs:
%	D - tesnor field (3x3xn1xn2xn3)
%	type -	1 (replace a non-spd tensor with the nearest spd tensor)
%			2 (force all the eigenvalues to be positive)
%
% Output:
%	Dnew - SPD tensor field (3x3xn1xn2xn3)
%
% Written by JK Gahm, UCLA. 01/17/2013.
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

if nargin<2
  type = 1;
end

RES = size(D);
RES = RES(3:end);
x = tenspd(D);
i = find(x==0);
[row,col,slc] = ind2sub(RES,i);
Dnew = D;

switch type
  case 1 % find the nearest spd tensor
    for j=1:length(i)
      nbs = 0;
      k = [];
      while isempty(k)
        nbs = nbs+1;
        x1 = x(max(1,row(j)-nbs):min(RES(1),row(j)+nbs),max(1,col(j)-nbs):min(RES(2),...
          col(j)+nbs),max(1,slc(j)-nbs):min(RES(3),slc(j)+nbs));
        D1 = D(:,:,max(1,row(j)-nbs):min(RES(1),row(j)+nbs),max(1,col(j)-nbs):min(RES(2),...
          col(j)+nbs),max(1,slc(j)-nbs):min(RES(3),slc(j)+nbs));
        k = find(x1);
      end
      l = k(randi(length(k)));
      RES1 = size(D1);
      D1 = reshape(D1,3,3,prod(RES1(3:end)));
      Dnew(:,:,row(j),col(j),slc(j)) = D1(:,:,l);
    end
  
  case 2 % use the absolute values of the eigenvalues
    [ev,ed] = teneig(D);
    for j=1:length(i)
      Dnew(:,:,i(j)) = ev(:,:,i(j))*abs(ed(:,:,i(j)))*ev(:,:,i(j))';
    end
end
      