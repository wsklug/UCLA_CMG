% This function calculates the tensor mode.
%
% SYNTAX:  M=ten_mode(D);
%          M=ten_mode(D,normdevD);
%
% INPUTS:  D - Tensor array with dimensions [3x3xLxMxN]
%          normdevD - [optional]  Tensor array with dimensions [LxMxN].
%          Sometimes normdevD is already computed in which case providing
%          it speeds up the computation.
%
% OUTPUTS: M - Tensor mode array with dimnesions [LxMxN]
%
% THEORY: M=det(D/ten_norm(D));
%
% DBE 2004.11.15
% DBE 2006.09.04 Updated to work with new ten_** functions.
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

function M=ten_mode(D,normdevD);

if nargin==1
  normdevD=ten_norm(ten_dev(D));
  ten_norm_dev=squeeze(repmat(permute(normdevD,[4 5 1 2 3]),[3 3 1 1 1]));
  M=3*sqrt(6)*ten_det(ten_dev(D)./ten_norm_dev);
elseif nargin==2
  ten_norm_dev=squeeze(repmat(permute(normdevD,[4 5 1 2 3]),[3 3 1 1 1]));
  %ten_norm_dev=repmat(permute(normdevD,[4 5 1 2 3]),[3 3 1 1 1]);
  M=3*sqrt(6)*ten_det(ten_dev(D)./ten_norm_dev);
else
  error('Too many input arguments.');
end
% warning('off');
%
% if nargin==1
%   ten_norm_dev=ten_norm(ten_dev(D));
%   ten_norm_dev=repmat(permute(ten_norm_dev,[4 5 1 2 3]),[3 3 1 1 1]);
%   ten_norm_dev=squeeze(ten_norm_dev);
%   M=3*sqrt(6)*ten_det(ten_dev(D)./ten_norm_dev);
% elseif nargin==2
%   if ndims(D)>2
%     M=zeros(sizes(D,3:ndims(D)));
%   else
%     M=0;
%   end
%   for j=1:size(normdevD,1)
%     for k=1:size(normdevD,2)
%       for l=1:size(normdevD,3)
%         M(j,k,l)=3*sqrt(6)*det(D(:,:,j,k,l)./normdevD(j,k,l));
%       end
%     end
%   end
% end
% warning('on');
%
% % M=3*sqrt(6)*det(D/ten_norm(D));

return