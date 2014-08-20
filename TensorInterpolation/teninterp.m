function Y = teninterp(D,X,interp_type)
%
% Diffusion tensor field interpolation
%
% Syntax: Y = teninterp(D,X,interp_type)
%
% Inputs:
%   D - tensor data (3x3xLxMxN)
%   X - interpolated (x,y,z) locations (nx3)
%   interp_type - tensor interpolation method
%                 'EU': Euclidean
%                 'LE': Log-Euclidean
%                 'Ki','Ri': Geodesic-loxodrome
%                 'NN': Nearest neighborhood
%                 'DY': Dyadic tensor for eigenvectors & linear invariant
%                   for eigenvalues
%
% Output:
%   Y - interpolated tensor component data (nx6)
%                 (D11,D12,D13,D22,D23,D33)
%
% Please cite the following when using this function:
%   Gahm JK et al. Linear invariant tensor interpolation applied to cardiac diffusion tensor MRI. MICCAI 2012, Part II. LNCS 7511, pp. 494-501, 2012.  
%   Gahm JK and Ennis DB. Dyadic tensor-based interpolation of tensor orientation: application to car- diac DT-MRI. STACOM, MICCAI Workshop 2013, LNCS 8330, pp. 135?142, 2014.
%
% Written by JK Gahm, UCLA. 02/06/2013.
% Ennis Lab @ UCLA; http://mrrl.ucla.edu


NUM_STEPS = 500;

n = size(X,1);
Y = zeros(n,6);
RES = size(D);
RES = RES(3:end);

% compute invariants for DY
if strcmp(interp_type,'DY')
  fprintf('Computing tensor invariants...\n')
  [K1,K2,K3,R1,R2,R3] = ten_ortho_invariants(D);
    
  fprintf('Spectral decomposition...\n')
  EV = teneig(D);
end

% compute matrix logarithm of tensor data
if strcmp(interp_type,'LE')
  T = reshape(D,[3 3 RES(1)*RES(2) RES(3)]);
  fprintf('Transforming tensor data into log-euclidean space')
  for i=1:RES(3)
    for j=1:RES(1)*RES(2)
      T(:,:,j,i) = logm(T(:,:,j,i));
    end
    fprintf('.')
  end
  fprintf('\n')
  D = reshape(T,[3 3 RES]);
  clear T
end

tic
for i=1:n
  
  % determine cubic end points for interpolation
  x1 = floor(X(i,:));
  x2 = x1+1;
  
  % tensor data at the end points
  D1 = D(:,:,x1(1):x2(1),x1(2):x2(2),x1(3):x2(3));
  D2 = ten2sym3(reshape(D1,[3 3 8])); % sym3 type
  
  % weighting function
  w = X(i,:)-x1;
  wf = [(1-w(1))*(1-w(2))*(1-w(3)) w(1)*(1-w(2))*(1-w(3)) ...
    (1-w(1))*w(2)*(1-w(3)) w(1)*w(2)*(1-w(3)) ...
    (1-w(1))*(1-w(2))*w(3) w(1)*(1-w(2))*w(3) ...
    (1-w(1))*w(2)*w(3) w(1)*w(2)*w(3)];
  
  switch interp_type
    case 'NN' % nearest neighbor
      w = X(i,:)-x1;
      id = zeros(3,1);
      % determine the nearest index for each dimension
      for j=1:3
        if w(j)<0.5
          id(j) = 1;
        elseif w(j)>0.5
          id(j) = 2;
        else
          id(j) = round(rand(1))+1; % otherwise chosen randomly
        end
      end
      Y(i,:) = ten2sym3(D1(:,:,id(1),id(2),id(3)));
    
    case 'DY' % dyadic tensor-based
      TR = K1(x1(1):x2(1),x1(2):x2(2),x1(3):x2(3));
      FA = R2(x1(1):x2(1),x1(2):x2(2),x1(3):x2(3));
      MD = R3(x1(1):x2(1),x1(2):x2(2),x1(3):x2(3));
      EV1 = EV(:,:,x1(1):x2(1),x1(2):x2(2),x1(3):x2(3));
      
      % tensor shape interpolation
      tr = sum(wf'.*TR(:));
      fa = sum(wf'.*FA(:));
      md = sum(wf'.*MD(:));
      ed = diag(ten_inv2lam([tr fa md],'combined')); 
      % tensor orientation interpolation
      ev = dyads_interp(EV1,w,ed);
      Y(i,:) = ten2sym3(ev*ed*ev');
      
    case 'EU' % Euclidean (linear)
      Y(i,:) = sum(repmat(wf,[6 1]).*D2,2);
  
    case 'LE' % log-Euclidean 
      Y(i,:) = ten2sym3(expm(sym32ten(sum(repmat(wf,[6 1]).*D2,2))));
      
    case {'Ki','Ri'} % geodesic-loxodrome
                     % *NOT available without teem library
      if strcmp(interp_type,'Ki')
        interp_type1 = 2;
      else
        interp_type1 = 3;
      end
      w1 = floor((X(i,:)-x1)*NUM_STEPS);
      Y(i,:) = tenInterpTeem(D2,NUM_STEPS,w1,interp_type1);  
  
    otherwise
      error('Wrong interpolation method!');  
  end  

  if mod(i,5000)==0
    fprintf('%1.0f%%\n',i/n*100)
  end
  
end

fprintf('total execution time = %.2f\n',toc)
