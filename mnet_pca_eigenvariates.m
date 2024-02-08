function [X,V,expvar, whitesig, dwM, whiteM] = mnet_pca_eigenvariates(Y,Nmodes,mode)
% [E,V,expvar] = mnet_pca_eigenvariates(Y,Nmodes,mode)
% Input:
%   Y: time x voxel
%   Nmodes: Number of components to decompose
%   mode:1 for spm_pca, 2 for other PCA
% Output:
%   X: eigenvalues
%   V: eigenvectors, voxel x num of eigenvalues
% by Hae-Jeong Park, Sep 20, 2014

if nargin<2, Nmodes=laplace_pca(Y); end
if nargin<3, mode=1; end
whitesig=[]; dwM=[];whiteM=[];
v=[]; expvar=[];
[m, n]   = size(Y); % Y=time x voxel

I=isfinite(Y); id=find(sum(I)~=0); 
if length(id)<1
    return;
else 
    Y=Y(:,id);
end
sz=size(Y,2);
I=isfinite(Y); idoff=find(sum(I')==sz); 
if length(idoff)<size(Y,1)*0.8, %too much infinite..
    Y(~I)=mean(Y(I));
    I=isfinite(Y); idoff=find(sum(I')==sz); 
end
Y=Y(idoff,:);

if isempty(Y), return; end
if mode(1)==1
    [X1,V,expvar]=spmsvd(Y,Nmodes); 
elseif mode(1)==2
    [whitesig, dwM, X1, V, whiteM, expvar] = spmsvd2(Y,Nmodes);
    X1=whitesig;
elseif mode(1)==3
    if numel(mode)==1
        norm=0;
    else
        norm=mode(2);
    end
    
    [V,X1,s] = pcasvd(Y,Nmodes,norm);s=diag(s); 
    V=V(1:Nmodes,:)';X1=X1(:,1:Nmodes);  
    expvar = cumsum(s)/sum(s);
end

X=nan(m,size(X1,2));
X(idoff,:)=X1;
end


function [L,V,expvar]=spmsvd(Y,Nmodes)
[m, n]   = size(Y); % Y=time x voxel
if nargin<1, Nmodes=m; end
if numel(Nmodes)==1,
    pcidx=1:Nmodes;
else pcidx=Nmodes; 
end

if m>n %for time series analysis in column, time size is vigger than voxel size
    [V, s, V] = svd(Y'*Y);
    s       = diag(s);
    V       = V(:,pcidx);
    u       = Y*V;    
    sq=sqrt(s(pcidx)); 
    for i=1:length(sq), u(:,i)=u(:,i)/sq(i); end
else
    [u, s, u] = svd(Y*Y');
    s       = diag(s);
    u       = u(:,pcidx);
    sq=sqrt(s(pcidx)); 
    V       = Y'*u;             
    for i=1:length(sq), V(:,i)=V(:,i)/sq(i); end
end  
d       = sign(sum(V));
sqn=sq/sqrt(n); 
L=u;
for i=1:length(sq), L(:,i)=u(:,i)*d(i)*sqn(i); end
expvar = cumsum(s)/sum(s); expvar=expvar'; % explained var; aka, condition numbers 
end

function [whitesig, dwM, L, V, whiteM, condnum] = spmsvd2(Y,Nmodes)
%function [yev,v,expvar]=spmsvd(Y,Nmodes)
[m, n]   = size(Y); % Y=time x voxel
if nargin<1, Nmodes=m; end
if numel(Nmodes)==1
    pcidx=1:Nmodes;
else pcidx=Nmodes; 
end

if m>n %for time series analysis in column, time size is vigger than voxel size
    %[v, s, v] = svd(Y'*Y);
    [v, s, v] = svd(cov(Y));
    %s       = diag(s);
    V=v;
    v       = v(:,pcidx);  
    sq=s(pcidx,pcidx); 
    
    whiteM=sqrtm(sq) \ v'; % use gaussian elimination
    dwM=v*sqrtm(sq);
    whitesig=Y * whiteM'; % PCs (=whitened signal)
    L=diag(s);
    %for i=1:length(sq), u(:,i)=u(:,i)/sq(i); end
else
    %[u, s, u] = svd(Y*Y');
    [v, s, v] = svd(cov(Y'));
    %s       = diag(s);
    V=v;
    v       = v(:,pcidx);
    sq=s(pcidx,pcidx); 
    
    whiteM=sqrtm(sq) \ v'; % use gaussian elimination
    dwM=v*sqrtm(sq);
    whitesig=Y' * whiteM'; % PCs (=whitened signal)  
    L=diag(s);
    %for i=1:length(sq), v(:,i)=v(:,i)/sq(i); end
end 
condnum = cumsum(L)/sum(L);

end
%d       = sign(sum(v));
%sqn=sq/sqrt(n); 
%yev=u;
%for i=1:length(sq), yev(:,i)=u(:,i)*d(i)*sqn(i); end
%end



% runpca -  perform principal component analysis (PCA) using singular value 
%           decomposition (SVD) using Matlab svd() or svds()
%                      >> inv(eigvec)*data = pc;
% Usage:
%             >> [pc,eigvec,sv] = runpca(data);
%             >> [pc,eigvec,sv] = runpca(data,num,norm)
%
%   data    input data matrix (rows are variables, columns observations)
%   num     number of principal comps to return  {def|0|[] -> rows in data}
%   norm    1/0 = do/don't normalize the eigvec's to be equivariant 
%                                                {def|0 -> no normalization}
% Outputs:
%   pc      the principal components, i.e.        >> inv(eigvec)*data = pc;
%   eigvec  the inverse weight matrix (=eigenvectors). >> data = eigvec*pc; 
%   sv      the singular values (=eigenvalues)

function [pc,M,S] = pcasvd(data,N,norm)

% Colin Humphries, CNL / Salk Institute, La Jolla CA
% 01/31/00 renamed runpca() and improved usage message -sm

BIG_N = 50; % for efficiency, switch to sdvs() when BIG_N<=N or N==rows

if nargin < 1
  help runpca
  return
end

rows = size(data,1);

if nargin < 3
  norm = 0;
elseif isempty(norm)
  norm = 0;
end

if nargin < 2
  N = 0;
end
if isempty(N)
  N = 0;
end

if N == 0  | N == rows
  N = rows;
  [U,S,V] = svd(data',0);   % performa SVD
  if norm == 0
    pc = U';
    M = (S*V')';
  else % norm
    pc = (U*S)';
    M = V;
  end
else
  if N > size(data,1)
    error('N must be <= the number of rows in data.')
  end
  if N <= BIG_N | N == rows
     [U,S,V] = svd(data',0);
  else
     [U,S,V] = svds(data',N);
  end
  if norm == 0
    pc = U';
    M = (S*V')';
  else % norm
    pc = (U*S)';
    M = V;
  end  
  if N > BIG_N & N < rows
    pc = pc(1:n,:);
    M = M(:,1:N);
  end
end
%S = diag(S(1:N,1:N));
end


function [k,p] = laplace_pca(data, e, d, n)
% Estimate latent dimensionality by Laplace approximation.
%
% k = LAPLACE_PCA([],e,d,n) returns an estimate of the latent dimensionality
% of a dataset with eigenvalues e, original dimensionality d, and size n.
% LAPLACE_PCA(data) computes (e,d,n) from the matrix data 
% (data points are rows)
% [k,p] = LAPLACE_PCA(...) also returns the log-probability of each 
% dimensionality, starting at 1.  k is the argmax of p.

%PMTKauthor Tom Minka
%PMTKurl http://research.microsoft.com/en-us/um/people/minka/papers/pca/

if ~isempty(data)
  [n,d] = size(data);
  m = mean(data);
  data0 = data - repmat(m, n, 1);
  e = svd(data0,0).^2;
end
e = e(:);
% break off the eigenvalues which are identically zero
i = find(e < eps);
e(i) = [];

% logediff(i) = sum_{j>i} log(e(i) - e(j))
logediff = zeros(1,length(e));
for i = 1:(length(e)-1)
  j = (i+1):length(e);
  logediff(i) = sum(log(e(i) - e(j))) + (d-length(e))*log(e(i));
end
cumsum_logediff = cumsum(logediff);

inve = 1./e;
% invediff(i,j) = log(inve(i) - inve(j))  (if i > j)
%               = 0                       (if i <= j)
invediff = repmat(inve,1,length(e)) - repmat(inve',length(e),1);
invediff(invediff <= 0) = 1;
invediff = log(invediff);
% cumsum_invediff(i,j) = sum_{t=(j+1):i} log(inve(t) - inve(j))
cumsum_invediff = cumsum(invediff,1);
% row_invediff(i) = sum_{j=1:(i-1)} sum_{t=(j+1):i} log(inve(t) - inve(j))
row_invediff = sum(cumsum_invediff);
% row_invediff(k) = sum_{i=1:(k-1)} sum_{j=(i+1):k} log(inve(j) - inve(i))

loge = log(e);
cumsum_loge = cumsum(loge);

cumsum_e = cumsum(e);

dn = length(e);
kmax = length(e)-1;
%dn = d;
%kmax = min([kmax 15]);
ks = 1:kmax;
% the normalizing constant for the prior (from James)
% sum(z(1:k)) is -log(p(U))
z = log(2) + (d-ks+1)/2*log(pi) - gammaln((d-ks+1)/2);
cumsum_z = cumsum(z);
for i = 1:length(ks)
  k = ks(i);
  %e1 = e(1:k);
  %e2 = e((k+1):length(e));
  %v = sum(e2)/(d-k);
  v = (cumsum_e(end) - cumsum_e(k))/(d-k);
  p(i) = -cumsum_loge(k) - (d-k)*log(v);
  p(i) = p(i)*n/2 - cumsum_z(k) - k/2*log(n);
  % compute h = logdet(A_Z)
  h = row_invediff(k) + cumsum_logediff(k);
  % lambda_hat(i)=1/v for i>k
  h = h + (d-k)*sum(log(1/v - inve(1:k)));
  m = d*k-k*(k+1)/2;
  h = h + m*log(n);
  p(i) = p(i) + (m+k)/2*log(2*pi) - h/2;
end
[pmax,i] = max(p);
k = ks(i);
  
%p(3)
%figure(1)
%plot(p)
end
