% This function constructs a tensor from a set of three invariants.
%
% SYNTAX: T=ten_inv2ten(K,TYPE)
%
% INPUTS:  K    - [1x3] or [3x1] vector of cylindrical (Ki) or spherical (Ri) invariants
%          type - 'cylindrical', 'Ki' -OR- 'Ri', 'spherical'
%
% OUTPUTS: T    - 3x3 resultant tensor
%
% DBE@UCLA 2010.03.11
% Ennis Lab @ UCLA; http://mrrl.ucla.edu
function lam=ten_inv2lam(K,type)

%% Define a numerical/floating point error bound
EPS=1e-6;

%% Define inputs if none are provided...
if nargin==0
  K(1)=1;    % NORM
  K(2)=0.5;  % FA
  K(3)=-1;   % MODE
  type='spherical';

%   K(1)=1.123;    % TRACE
%   K(2)=0.67;  % NORM_DEV
%   K(3)=-sqrt(2)/2;   % MODE
%   type='cylindrical';
end

%% Perform a check on the input invariant values...
lam = zeros(size(K));

switch lower(type)
  case {'spherical','ri'}
    i = K(:,1)<0-EPS | K(:,2)<0-EPS | K(:,2)>1+EPS | K(:,3)<-1-EPS | K(:,3)>1+EPS;
    K(i,1) = NaN; K(i,2) = NaN; K(i,3) = NaN;
    theta = -(1/3)*asin(K(:,3));
    lam(:,1) = (1/3)*K(:,1).*sqrt(3-2*K(:,2).^2)+(2/3)*K(:,1).*K(:,2).*sin(theta+(2/3)*pi);
    lam(:,2) = (1/3)*K(:,1).*sqrt(3-2*K(:,2).^2)+(2/3)*K(:,1).*K(:,2).*sin(theta         );
    lam(:,3) = (1/3)*K(:,1).*sqrt(3-2*K(:,2).^2)+(2/3)*K(:,1).*K(:,2).*sin(theta-(2/3)*pi);
  
  case {'cylindrical','ki'}
    i = K(:,1)<0-EPS | K(:,2)<0-EPS | K(:,3)<-1-EPS | K(:,3)>1+EPS;
    K(i,1) = NaN; K(i,2) = NaN; K(i,3) = NaN;
    theta = -(1/3)*asin(K(:,3));
    lam(:,1) = (1/3)*K(:,1)+sqrt(2/3)*K(:,2).*sin(theta+(2/3)*pi);
    lam(:,2) = (1/3)*K(:,1)+sqrt(2/3)*K(:,2).*sin(theta         );
    lam(:,3) = (1/3)*K(:,1)+sqrt(2/3)*K(:,2).*sin(theta-(2/3)*pi);
    
  case 'combined' % K1,R2,R3
    i = K(:,1)<0-EPS | K(:,2)<0-EPS | K(:,2)>1+EPS | K(:,3)<-1-EPS | K(:,3)>1+EPS;
    K(i,1) = NaN; K(i,2) = NaN; K(i,3) = NaN;
    theta = (1/3)*acos(K(:,3));
    lam(:,1) = (1/3)*K(:,1)+(2/3)*K(:,1).*K(:,2).*cos(theta)./sqrt(3-2*K(:,2).^2);
    lam(:,2) = (1/3)*K(:,1)+(2/3)*K(:,1).*K(:,2).*cos(theta-(2/3)*pi)./sqrt(3-2*K(:,2).^2);
    lam(:,3) = (1/3)*K(:,1)+(2/3)*K(:,1).*K(:,2).*cos(theta+(2/3)*pi)./sqrt(3-2*K(:,2).^2);    
  
  otherwise
    error('TYPE not understood.');
end
