function mnet_savedcms(DCMs)
for i=1:length(DCMs)
    try DCM=DCMs{i};
    catch, DCM=DCMs(i); 
    end
    if ~isempty(DCM.name)
        save(DCM.name,'DCM',spm_get_defaults('mat.format'));
    end
end
end