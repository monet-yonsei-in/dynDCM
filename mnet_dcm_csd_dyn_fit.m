function [CSD,PEB0,M0]=mnet_dcm_csd_dyn_fit(y,TR,roinames,wndsz,options,dcmfile)
% [CSD,PEB0,M0]=mnet_dcm_csd_dyn_fit(y,TR,roinames,wndsz,options,dcmfile)
% Convert DCM diagnoal terms to be real Hz
% INPUT:
% OUTPUT:
% Hae-Jeong Park @ Yonsei University, Aug. 15, 2017
if nargin<2, TR=2; end
if nargin<3, roinames=[]; end
if nargin<4, wndsz=[]; end %[window size overlap size]
if nargin<5, options=[]; end
if nargin<6, dcmfile=[]; end
if isempty(wndsz), wndsz=128; end
if ~isfield(options,'nopeb'), options.nopeb=1; end
if ~isfield(options,'peb_niter'), options.peb_niter=4; end
if ~isfield(options,'peb_skipfirst'), options.peb_skipfirst=0; end
if ~isfield(options,'delete_temp_files'), options.delete_temp_files=1; end
if ~isfield(options,'saveall'), options.saveall=1; end
if ~isfield(options,'parallel'), options.parallel=true; end
if ~isfield(options,'hanning'), options.hanning=0; end

DCM={}; PEB0=[];M0=[];
if isnumeric(y)    
    y=double(y);
    dcm=mnet_dcm_csd_model(y,TR,roinames,options,dcmfile);

    if length(wndsz)==1, wndsz=[wndsz 0]; end
    nscan=size(y,1);
    try
        offset=wndsz(1)-wndsz(2);
        nwnd=floor((nscan-wndsz(2))/offset);        
    catch
        return; 
    end
    if options.hanning>0
        hw=spm_hanning(wndsz(1));
    end
    for s = 1:nwnd
        if options.hanning>0
            dcm.Y.y=hw.*y([1:wndsz(1)]+(s-1)*offset,:);   
        else 
            dcm.Y.y=y([1:wndsz(1)]+(s-1)*offset,:);   
        end
        [p,f,e]=fileparts(dcmfile);
        dcmfn=fullfile(p,[f '_' num2str(s) '.mat']);
        dcm.name=dcmfn;    
        DCM{s,1}   = dcm;
    end
else
    DCM=y;
    nwnd=length(DCM);
    for i=1:nwnd
        if ~exist(DCM{i,1}.name,'file')
            mnet_savedcms(DCM); break;
        end
    end
end

if options.nopeb
    CSD  = spm_dcm_fit(DCM,options.parallel); M0=[]; PEB0=[];
else
    [CSD,PEB0,M0] = spm_dcm_peb_fit_par(DCM,[],[],options.peb_niter,options.peb_skipfirst);
end

if options.delete_temp_files
    for s = 1:nwnd
        dcmfn=DCM{s,1}.name;
        delete(dcmfn);
    end
end

[p,f,e]=fileparts(dcmfile); 
dcmfile1=fullfile(p,[f '_dyn.mat']);    
save(dcmfile1,'CSD','M0','PEB0','roinames','wndsz','options','nscan','nwnd','-v7.3');    

end
