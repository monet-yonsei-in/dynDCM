function [M,V]=mnet_dcm_estimate_dynamics(CSD,options)
% [M,V]=mnet_dcm_estimate_dynamics(CSD,options)
% Convert DCM diagnoal terms to be real Hz
% INPUT:
% OUTPUT:
% Hae-Jeong Park @ Yonsei University, Aug. 15, 2017

if nargin<2, options.mode=0; end
if isnumeric(options), opt=options; options=[]; 
    if ~isempty(opt),options.mode=opt(1); options.ndct=opt(2); end; 
end
try options.display, catch options.display=0; end
try beta=options.beta; catch, beta  = 32; end
try alpha=options.alpha; catch, alpha  = 16; end
try eta=options.eta; catch, eta  = 0; end
try options.dcm2hz; catch options.dcm2hz=0; end
V=[];
nwnd=length(CSD);
if options.mode==0 %AVG
    X=ones(nwnd,1);
elseif options.mode==1 %DCT
    mtx=spm_dctmtx(nwnd);    
    if numel(options.ndct)==1,
        X=[ones(nwnd,1) mtx(:,2:options.ndct)]; 
    else
        id=find(options.ndct==1); options.ndct(id)=[];
        X=[ones(nwnd,1) mtx(:,options.ndct)]; 
    end
    
elseif options.mode==2 % PCA
    Y=zeros(length(CSD),prod(size(CSD{1,1}.Ep.A))); mY=[];
    for i=1:length(CSD)
       a=CSD{i,1}.Ep.A;
       if options.dcm2hz, a=mnet_dcm2hz(a); end
       Y(i,:)=a(:)'; 
       mY=mean(Y);
    end
    [mtx,V] = mnet_pca_eigenvariates(Y,options.ndct,1);
    if ~isempty(mY)
        V(:,2)=mY'; V=V(:,[2 1]);
    end
    X=[ones(nwnd,1) mtx];
 elseif options.mode==3 % supfPCA
    Y=zeros(length(CSD),prod(size(CSD{1,1}.Ep.A)));
    for i=1:length(CSD)
       a=CSD{i,1}.Ep.A;
       if options.dcm2hz, a=mnet_dcm2hz(a); end
       Y(i,:)=a(:)'; 
    end
    try r=options.ndct(2); catch r=4; end
    mY=mean(Y);Y=Y-repmat(mY,nwnd,1); %npca=npca-1;
    [mtx,b,avgx] = mnet_supfpca(Y',options.ndct(1),r); %ones(length(CSD),1)
    X=mtx; %[avgx' mtx]; 
    X=[ones(length(CSD),1) spm_detrend(X)];
elseif options.mode==4 % PCA
    Y=zeros(length(CSD),prod(size(CSD{1,1}.Ep.A)));
    for i=1:length(CSD)
       a=CSD{i,1}.Ep.A;
       if options.dcm2hz, a=mnet_dcm2hz(a); end
       Y(i,:)=a(:)'; 
    end
    demean=1;
    npca=options.ndct(1);
    if numel(options.ndct)>1
        demean=options.ndct(2);
    end
    mY=[];
    if demean
        mY=mean(Y);Y=Y-repmat(mY,nwnd,1); npca=npca-1;
    end
    [mtx,V] = mnet_pca_eigenvariates(Y,npca,2);
    if ~isempty(mY)
        V(:,2)=mY'; V=V(:,[2 1]);
    end
    X=[ones(nwnd,1) mtx];      
elseif options.mode==5
    X=options.X;    
else
    X=ones(nwnd,1);
end

if options.display
    figure;h=plot(X);
    set(gca,'YLim',[-1.5 1.5],'linewidth',1.5,'FontSize',14); set(h,'Linewidth',2);set(gcf,'color',[1 1 1]);
    title('Estimated Dynamics'); xlabel('Window');
end


M.X   = X;                  % between session explanatory variables
M.hE  = eta;                         % prior expectation of log precision
M.hC  = 1/alpha;                   % prior covariance of precision
M.bE  = CSD{1}.M.pE;               % prior expectations over sessions
M.bC  = CSD{1}.M.pC;               % prior covariance over sessions
M.pC  = CSD{1}.M.pC/beta;          % prior covariance between sessions
end