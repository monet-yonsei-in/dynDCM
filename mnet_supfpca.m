function [V_sf1,B_sf1,avgX]=mnet_supfpca(X,Y,r)
% X=[channel x time]
if nargin<2, Y=[]; end
if nargin<3, r=[]; end

V_sf1=[];B_sf1=[];

sz=size(X);
avgX=mean(X,1); %effective connectivity matrix
%X=X-repmat(avgX,sz(1),1);
%X=X-mean(avgX(:));
if isempty(r)
    [~,D,~]=svd(X,'econ');
    scree=diag(D).^2/sum(diag(D).^2);
    figure;
    plot([1:10],scree(1:10),'b.-','linewidth',1.5,'markersize',30);
    xlabel('Rank','fontsize',15);
    ylabel('Variance Explained','fontsize',15);
    set(gca,'fontsize',15);
    title('Scree Plot of X','fontsize',20)
    r=input('Input r:\n');
end


ndct=floor(sz(1)*2/3);
if numel(Y)==1, ndct=Y; Y=[]; end
if isempty(Y)
    Y=spm_dctmtx(sz(1),ndct);
    %Y=Y(:,2:end);
    %Y=Y-repmat(mean(Y,1),sz(1),1); % center supervision data
end

[B_sf1,V_sf1,U_sf1,se2_sf1,Sf_sf1]=SupSFPCA(Y,X,r,struct('lambda',0,'gamma',0));

[~,~,V]=svds(X,r);
V=bsxfun(@times,V,[-1,ones(1,r-1)]);

signV=diag(sign(V_sf1'*V))';
V_sf1=bsxfun(@times,V_sf1,signV);
B_sf1=bsxfun(@times,B_sf1,signV);

end


function [B,V,U,se2,Sf]=SupSFPCA(Y,X,r,paramstruct)
% This function conducts supervised sparse and functional principal
% component analysis by fitting the SupSVD model
%       X=UV' + E 
%       U=YB + F  
% where X is an observed primary data matrix (to be decomposed), U is a latent score
% matrix, V is a loading matrix, E is measurement noise, Y is an observed
% auxiliary supervision matrix, B is a coefficient matrix, and F is a
% random effect matrix.
%
% It decomposes the primary data matrix X into low-rank
% components, while taking into account many different features: 1)
% potential supervision from any auxiliary data Y measured on the same
% samples; 2) potential smoothness for loading vectors V (for functional
% data); 3) sparsity in supervision coefficients B and loadings V (for variable
% selection).
%
% It is a very general dimension reduction method that subsumes
% PCA, sparse PCA, functional PCA, supervised PCA, etc as special cases.
% See more details in 2016 JCGS paper "Supervised sparse and 
% functional principal component analysis" by Gen Li, Haipeng Shen, and
% Jianhua Z. Huang.
%
% Input:
%   Y       n*q (column centered) auxiliary data matrix, rows are samples and columns are variables
%
%   X       n*p (column centered) primary data matrix, which we want to decompose. 
%           rows are samples (matched with Y) and columns are variables
%
%   r       positive scalar, prespecified rank (r should be smaller than n and p)
%
%   paramstruct
%       lambda      0 or 1 (default=1, sparse loading), sparsity index for loadings
%     
%       alpha       0 or 1 (default=1, smooth loading), smoothness index for loadings
%
%       gamma       0 or 1 (default=1, sparse coefficient), sparsity index
%                   for supervision coefficients. Note: if gamma is set to be 0,
%                   Y must have q<n to avoid overfitting; if gamma is set to be 1,
%                   then it can handle high dimensional supervision Y
%
%       Omega       p*p symmetric positive semi-definite matrix for
%                   smoothness penalty (default is for evenly spaced data)
%                   Note: only change this if you have unevenly spaced
%                   functional data X
%
%       convg_thres    positive scalar (default=1E-6), overall convergence threshold
%
%       vconvg_thres   positive scalar (default=1E-4), convergence
%                      threshold for the proximal gradient descent algorithm  
%                      for estimating each loading vector 
%
%       max_niter   scalar (default=1E3), max number of overall iterations
%
%       vmax_niter  scalar (default=1E2), max number of iterations for estimating
%                   each loading vector 
%
% Output:
%   B       q*r coefficient matrix of Y on the scores of X, 
%           maybe sparse if gamma=1
%
%   V       p*r loading matrix of X, each column has norm 1, but no strict orthogonality
%           because of sparsity and smoothness. If lambda=1, V is sparse;
%           if alpha=1, each column of V is smooth
%
%   U       n*r score matrix of X, conditional expectation of random scores, 
%           no strict orthogonality
%
%   se2     scalar, variance of measurement error in the primary data X
%
%   Sf      r*r diagonal covariance matrix, for random effects (see paper)
%
% Note: Essentially, U and V are the most important output for dimension
% reduction purpose as in PCA or SVD. 
%
% %%%%%%%%%%%%%%%%%%%%
% Cite: 
% @article{li2015supervised,
%   title={Supervised sparse and functional principal component analysis},
%   author={Li, Gen and Shen, Haipeng and Huang, Jianhua Z},
%   journal={Journal of Computational and Graphical Statistics},
%   number={just-accepted},
%   pages={00},
%   year={2015},
%   publisher={Taylor \& Francis}
% }
%
% Contact: Gen Li, PhD
%          Assistant Professor of Biostatistics, Columbia University
%          Email: gl2521@columbia.edu  
%
% CopyRight all reserved
% Last updated: 4/16/2016
% %%%%%%%%%%%%%%%%%%%%%%

if max(abs(mean(X,1)))>0.01 || max(abs(mean(Y,1)))>0.01
    %error('Columns of X and Y are not centered. exit...');
end;

ind_lam=1; % default: v sparse
ind_alp=1; % default: v smooth
ind_gam=1; % default: b sparse
ind_Omg=1; % default: Omega even spaced
max_niter=1e3; % max num of iter for EM
convg_thres=1E-6;  % EM convergence rule
vmax_niter=1e2; % max num of iter for each v_k
vconvg_thres=1E-4; % v_k convergence rule

if nargin > 3 ;   %  then paramstruct is an argument
  if isfield(paramstruct,'lambda') ;    %  then change to input value
    ind_lam = getfield(paramstruct,'lambda') ; 
  end ;
  if isfield(paramstruct,'alpha') ;    %  then change to input value
    ind_alp = getfield(paramstruct,'alpha') ; 
  end ;
  if isfield(paramstruct,'gamma') ;    %  then change to input value
    ind_gam = getfield(paramstruct,'gamma') ; 
  end ;
  if isfield(paramstruct,'Omega') ;    %  then change to input value
    Omega = getfield(paramstruct,'Omega') ; 
    ind_Omg=0;
  end ;
  if isfield(paramstruct,'convg_thres') ;    %  then change to input value
    convg_thres = getfield(paramstruct,'convg_thres') ; 
  end ; 
  if isfield(paramstruct,'vconvg_thres') ;    %  then change to input value
    vconvg_thres = getfield(paramstruct,'vconvg_thres') ; 
  end ; 

end;


[n,p]=size(X);
[n1,q]=size(Y);
if(rank(Y)<q)
    error('Do not run this code! Change initial and BIC df....');
    error('gamma cannot be set to zero!');
end;

% Pre-Check
if (n~=n1)
    error('X does not match Y! exit...');
elseif (rank(Y)~=q)
    error('Columns of Y are linearly dependent! exit...');
elseif (r>n || r>p)
    error('Too greedy on ranks! exit...');
end;

% set Omega
if(ind_Omg==1)
    Q=eye(p)*(-2);
    Q=spdiags(ones(p,1),1,Q);
    Q=spdiags(ones(p,1),-1,Q);
    Q=Q(:,2:(end-1));
    R=eye(p-2)*(2/3);
    R=spdiags(ones(p-2,1)*(1/6),1,R);
    R=spdiags(ones(p-2,1)*(1/6),-1,R);
    Omega=Q*inv(R)*Q'; % p*p
end;
oeig=eigs(Omega,1); % largest eig value of Omega


% initial est
[U,D,V]=svds(X,r); % initial V
U=U*D; % initial U
E=X-U*V';
se2=var(E(:)); % initial se2
B=inv(Y'*Y)*Y'*U; % initial B
Sf=diag(diag((1/n)*(U-Y*B)'*(U-Y*B))); % initial Sf
clear E D;




diff=1; % EM criterion
niter=0; % number of EM iter
while (niter<=max_niter && diff>convg_thres )
    % record last iter
    se2_old=se2;
    Sf_old=Sf;
    V_old=V;
    B_old=B;
    
    
    % E step
    % some critical values
    Sfinv=inv(Sf);
    weight=inv(eye(r)+se2*Sfinv); % r*r   
    cond_Mean=(se2*Y*B*Sfinv + X*V)*weight; % E(U|X), n*r
    cond_Var=Sf*(eye(r)-weight); % cov(U(i)|X), r*r   
    cond_quad=n*cond_Var + cond_Mean'*cond_Mean; % E(U'U|X), r*r
     
    
    % M step
    % estimate B and Sf
    if(ind_gam~=0)
        for k=1:r
            % Attention: lasso fcn always center X, Y, and return a
            % separate column of intercept; by default, it standardize each
            % column of X
            % Therefore, when using lasso, we always center the columns of
            % X and Y to avoid the trouble of intercept
            [SpParam,FitInfo]=lasso(Y,cond_Mean(:,k),'LambdaRatio',0,'Standardize',false); 
            BIC_score=n*log(FitInfo.MSE)+log(n)*FitInfo.DF;
            [~,ind]=min(BIC_score);
%             figure(1);clf;plot(FitInfo.Lambda,BIC_score);
%             figure(2);clf;plot(FitInfo.Lambda,SpParam);
            B(:,k)=SpParam(:,ind);
        end;
    else % if gamma=0, no B-sparsity
        B=inv(Y'*Y)*Y'*cond_Mean;
    end;
    %
    % estimate Sf
    Sf=diag(diag( (cond_quad + (Y*B)'*(Y*B)- (Y*B)'*cond_Mean- cond_Mean'*(Y*B) )/n )); % r*r
    %
    % estimate V
    for k=1:r % kth layer
        % some critical values
        theta=X'*cond_Mean(:,k)-(V_old*cond_quad(:,k)-V_old(:,k)*cond_quad(k,k)); % p*1
        c=cond_quad(k,k); % E(u_k'u_k|X), 1*1
        
        % select smooth tuning (LOOCV w/o refitting)
        if(ind_alp~=0)
            alphavec=0:0.1:10; % smooth tuning range
            cv_score=zeros(size(alphavec));
            for ialp=1:length(alphavec)
                alpha=alphavec(ialp);
                hat=inv(eye(p)+alpha*Omega);
                vest=hat*theta/c;
                cv_score(ialp)=(1/p)*sum(((theta/c-vest)./(1-diag(hat))).^2);
            end;
            [~,I]=min(cv_score);
            optalp=alphavec(I); % optimal alpha for this iteration
%           figure(k);clf;plot(alphavec,cv_score,'.-');title('Smoothness Tuning');
        else % no V-smoothness
            optalp=0;
        end;
 
        % specify sparsity tuning (for gaussian error)
        if(ind_lam~=0)
            optlam=sqrt(2*log(p)*se2_old/c);               
        else % no V sparsity
            optlam=0;
        end;
        
        
        L=1+optalp*oeig; 
        vk_old=V_old(:,k); % initial value for v_k is from last EM iteration
        vdiff=1;
        vniter=0;
        while(vniter<=vmax_niter && vdiff>vconvg_thres ) % iteration for estimating v_k
            df= -theta/c + (eye(p)+optalp*Omega)*vk_old;
            vk=soft_thr(vk_old-(1/L)*df, optlam/L);
            % set norm=1
            if norm(vk)==0
                warning('zero vector v!');
            else
                vk=vk/norm(vk);
            end;
            vdiff=norm(vk-vk_old)^2;
            vk_old=vk;
            vniter=vniter+1;
        end;
        V(:,k)=vk;
%         figure(k);clf;plot(vk);
    end;
    %         
    % Estimate se2
    se2=(trace(X*(X'-2*V*cond_Mean')) + n*trace(V'*V*cond_Var) + trace((cond_Mean'*cond_Mean)*(V'*V)))/(n*p);

        
    % stopping rule
    diff=norm(V-V_old,'fro')^2;
    niter=niter+1;

end;


% Print convergence information
if niter<max_niter
    disp(['SupFPCA converges at precision ',num2str(convg_thres),' after ',num2str(niter),' iterations.']);
else
    disp(['SupFPCA NOT converge at precision ',num2str(convg_thres),' after ',num2str(max_niter),' iterations!!!']);
end;



% reorder V and others
[~,I]=sort(diag(V'*X'*X*V),'descend');
V=V(:,I);
B=B(:,I);
Sf=Sf(I,I);


% output U
Sfinv=inv(Sf);
weight=inv(eye(r)+se2*Sfinv); % r*r
U=(se2*Y*B*Sfinv + X*V)*weight;

end

function out = soft_thr(in,lambda)
    out = sign(in).*max(abs(in) - lambda,0);
end
