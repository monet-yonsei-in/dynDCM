function DCM=mnet_dcm_csd_model(y,TR,roinames,options,dcmfile) 
% DCM=mnet_dcm_csd_model(y,TR,roinames,options,dcmfile)
% y: time x region
% TR: time of repetition (sampling time for scan)
% roinames: {'ROI1','ROI2',...};
% options.A: logical connection model, decide which parameter will be included in the estimation, the parameters of nonzero A, defaults: [] 
% options.pC: a priori for the covariance, constraining precision of the connection matrix, defaults:[] 
% options.A0 : initial point of the A matrix in the optimization
% options.maxnodes : if you have more than 8 regions then change this to
% suit your number of regions [default: 8]
% options.Fdcm : frequency ranges default=[0.0078 min(0.2, 1/(2*TR))] Hz
% options.globalfluctuation : if 0, each node flucutate, otherwise all nodes share a same parameter
% options.ARols :  if 1, conduct AR estimation using OLS otherwise using %Bayesian approach
% options.nfreqbins; number of frequency bins [default: nfreqbins=64]
% 
% dcmfile for DCM
% by Hae-Jeong Park, Ph.D. at May 13, 2016
% -------------------------------------------------------------------------

if nargin<2, TR=2; end
if nargin<3, roinames=[]; end
if nargin<4, options=[]; end
if nargin<5, dcmfile=[]; end

[T, n] = size(y); % T = number of oberservations; n =number of regions

% options
% -------------------------------------------------------------------------
options.induced     = 1; %important
options.two_state  = 0; %important in configuration of priors

%-- below is not so important... use default
options.nonlinear  = 0;
options.stochastic = 0;
options.nograph    = 1;
options.analysis   = 'CSD';    


try options.Fdcm; catch, options.Fdcm=[0.0078 min(0.2, 1/(2*TR))]; end
try options.globalfluctuation; catch, options.globalfluctuation=1; end
try options.ARols; catch, options.ARols=0; end
try options.nfreqbins; catch, options.nfreqbins=32; end
try options.maxnodes; catch, options.maxnodes=8; end


% response
% -------------------------------------------------------------------------
%DCM.Y.X0 %nuisances...
DCM.Y.y  = y;
DCM.Y.dt = TR;

if isempty(roinames)
    for i=1:n
        roinames{i}=sprintf('ROI%d',i);
    end
end  

DCM.Y.name=roinames;
% experimental inputs 
% -------------------------------------------------------------------------

DCM.U.u    = zeros(T,1);
DCM.U.name = {'null'};
DCM.U.dt=TR;


Ai=[];initA=[];pCov=[];
if isfield(options,'A'), Ai=options.A; options=rmfield(options,'A'); end
if isfield(options,'A'), initA=options.A0; options=rmfield(options,'A0'); end
if isfield(options,'A'), pCov=options.pC; options=rmfield(options,'pC'); end

if isempty(Ai),Ai   = ones(n,n); else Ai = logical(Ai); end
B   = zeros(n,n,0); 
C   = zeros(n,1);
D   = zeros(n,n,0);

[pE,pC]  = mnet_dcm_fmri_priors(Ai,B,C,D,options,pCov);
 
% this setup will be redone in spm_dcm_fmri_csd...
DCM.a    = logical(Ai);
DCM.b    = zeros(n,n,0);
DCM.c      = zeros(n,1);
DCM.d    = zeros(n,n,0);
DCM.options = options;

%% check for pre-specified priors
%--------------------------------------------------------------------------
if ~isempty(pCov)
    DCM.M.pE = pE;  DCM.M.pC = pC; 
end
% -------------------------------------------------------------------------
%% assign initial values for optimization
if ~isempty(initA)
    DCM.options.P = initA;
end
%%

if isempty(dcmfile),
    DCM.name = sprintf('DCM_%s_%s',date,mnet_time);
else
    DCM.name=dcmfile; 
    save(DCM.name,'DCM');
end


function str=mnet_time
a=fix(clock);str=[int2str(a(2)) int2str(a(3)) int2str(a(4)) int2str(a(5)) int2str(a(6))];

