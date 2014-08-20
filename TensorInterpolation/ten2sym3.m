function y = ten2sym3(x)
% Convert 3x3 symmetric tensors into 6-tensor component vectors
%
% Syntax: y = ten2sym3(x)

% Input:
%	x - 3x3 symmetric tensors (3x3xn1xn2xn3)
%
% Output:
%	y - 6-tensor component vectors (6xn1xn2xn3)
%
% Written by JK Gahm, UCLA. 01/17/2013.
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

n1 = size(x,3);
n2 = size(x,4);
n3 = size(x,5);

x = reshape(x,[9 n1*n2*n3]);
y = zeros(6,n1,n2,n3);
for i=1:n1*n2*n3
    %y(:,i) = [x(1,i);x(4,i);x(7,i);x(5,i);x(8,i);x(9,i)];
    y(:,i) = [x(1,i);x(2,i);x(3,i);x(5,i);x(6,i);x(9,i)];
end
