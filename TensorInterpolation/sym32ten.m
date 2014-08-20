function y = sym32ten(x)
% Convert 6-tensor component vectors into 3x3 symmetric tensors
%
% Syntax: y = sym32ten(x)

% Input:
%	x - 6-tensor component vectors (6xn)
%
% Output:
%	y - 3x3 symmetric tensors (3x3xn)
%
% Written by JK Gahm, UCLA. 01/17/2013.
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

n = size(x);

y = zeros([3 3 n(2:end)]);
for i=1:prod(n(2:end))
  if n(1)==6
    y(:,:,i) = [x(1,i) x(2,i) x(3,i); x(2,i) x(4,i) x(5,i); x(3,i) x(5,i) x(6,i)];
  elseif n(1)==7
    y(:,:,i) = [x(2,i) x(3,i) x(4,i); x(3,i) x(5,i) x(6,i); x(4,i) x(6,i) x(7,i)];
  end
end
