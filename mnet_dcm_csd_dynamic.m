function [PEB,BMA,CSD,PEB0,DCMa]=mnet_dcm_csd_dynamic(y,TR,roinames,wndsz,options,dcmfile,modeltitle)
% [PEB,BMA,CSD,PEB0,DCMa]=mnet_dcm_csd_dynamic(y,TR,roinames,wndsz,options,dcmfile,initialmode)
% INPUT:
% OUTPUT:
% Hae-Jeong Park @ Yonsei University, Aug. 15, 2017

if nargin<2, TR=2; end
if nargin<3, roinames=[]; end
if nargin<4, wndsz=[]; end %[window size overlap size]
if nargin<5, options=[]; end
if nargin<6, dcmfile=[]; end
if nargin<7, modeltitle=[]; end

PEB=[];BMA=[];PEB0=[];M=[]; CSD=[];DCMa=[];
initialestimate=1; skipmode=1;
if isnumeric(options) 
    opt=options; options=[]; 
    if ~isempty(opt)
        if numel(opt)<3
            options.mode=opt(1); 
            options.ndct=opt(2);
        else
            options.X=opt;
            options.mode=5;
            if isempty(modeltitle)
                options.ndct='hmmdv';
            else
                options.ndct='';
            end
        end
    end
end
try field=options.fields; catch, field = {'A'}; end %{'A','a'}                 % parameters of interest
try options.skipmode, skipmode=options.skipmode; catch, skipmode=1; end
try options.initialmode; initialmode=options.initialmode; catch, initialmode=1; end
try options.nopeb;  catch, options.nopeb=0; end
try options.bmc_peb; catch, options.bmc_peb = 0; end
try options.peb_bmc; catch, options.peb_bmc = 0; end
try options.total_epoch; catch, options.total_epoch=0; end

if initialmode>1, initialestimate=0; end %only for improving parallel computing resources

if ~isempty(dcmfile)
    [p,f,e]=fileparts(dcmfile); 
    dcmfile1=fullfile(p,[f '_dyn.mat']);    
    if exist(dcmfile1,'file') && skipmode
        initialestimate=0;
    else 
        initialestimate=1;
    end
else 
    dcmfile=datestr(now,30); id=strfind(dcmfile,':'); dcmfile(id)=[]; id=dcmfile==' '; dcmfile(id)='-';
    dcmfile=fullfile(pwd,['DCM_' deblank(dcmfile) '_' num2str(randi(1000,1)) '.mat']);
    [p,f,e]=fileparts(dcmfile); 
    dcmfile1=fullfile(p,[f '_dyn.mat']);
end

if ischar(y)
   for i=1:size(y,1)
      dcm=load(deblank(y(i,:)));
      CSD{i,1}=dcm.DCM;
   end
   initialestimate=0;
end
if iscell(y)
    CSD=y; initialestimate=0;
end

if initialestimate
    [CSD,PEB0,M0]=mnet_dcm_csd_dyn_fit(y,TR,roinames,wndsz,options,dcmfile);
    save(dcmfile1,'CSD','M0','PEB0','field','roinames','wndsz','-v7.3');
    if initialmode==1
        return;
    end    
else
    if ~isempty(dcmfile1)
        if exist(dcmfile1,'file'), load(dcmfile1); end
    end
end
nwnd=length(CSD);
try options.mode; catch, options.mode=1; end
try options.ndct; catch, options.ndct=3; end

modestr='';
[p,f,e]=fileparts(dcmfile); 

if options.mode==1 %DCT
    modestr=sprintf('DCT%s%dof%d',modeltitle,options.ndct,nwnd);
elseif options.mode==2 % PCA spm
    modestr=sprintf('PCA%s%dof%d',modeltitle,options.ndct,nwnd);
elseif options.mode==3 % supfPCA
    modestr=sprintf('FPCA%s%d_%d_W%d',modeltitle,options.ndct(2),options.ndct(1),nwnd);
elseif options.mode==4 %PCA flush
    modestr=sprintf('PCAFLUSH%s%dW%d',modeltitle,options.ndct(1),nwnd);
elseif options.mode==5
    modestr=sprintf('%s%sW%d',modeltitle,options.ndct,nwnd);
else
end

dcmfile2=fullfile(p,[f '_dyn_' modestr '.mat']);
if ~exist(dcmfile2,'file') 
    
    M=mnet_dcm_estimate_dynamics(CSD,options);
    for k=1:3
        PEB   = spm_dcm_peb(CSD,M,field);  % empirical Bayesian inversion
        if ~isempty(PEB) 
            if isnan(PEB.F)==0
                break;
            end
        end
        M.bC=M.bC*2;
    end
    if ~isempty(dcmfile2)
        save(dcmfile2,'PEB','M','options','field','-v7.3');
    end

    if options.bmc_peb %%  reduce in the 2nd levels...
        [BMC,PEBc]= spm_dcm_bmc_peb(CSD,M,'A'); 
        if exist(dcmfile2,'file')
            save(dcmfile2,'BMC','PEBc','-append');
        end 
    end
    if options.peb_bmc 
        BMA   = spm_dcm_peb_bmc(PEB);  
        if exist(dcmfile2,'file')
            save(dcmfile2,'BMA','-append');
        end 
    else BMA=[];
    end    % Bayesian model reduction and averaging
    
end


if options.total_epoch
    if exist(dcmfile,'file')
        load(dcmfile);
        try, DCM.Ep.A; catch
            DCMa=spm_dcm_fmri_csd(dcmfile);
        end
        delete(dcmfile);
        save(dcmfile1,'DCMa','-append');
    end    
end