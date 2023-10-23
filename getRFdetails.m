
% 110221: modified to rfStatsPix=[]. Not deleted to keep compatibility.

function [rfStatsDeg,rfStatsPix,LFPElectrodeList,EcogElectrodeList,subjectNames] = getRFdetails(subjectNamesA,folderName)

filename  = fullfile(folderName,'RFDetails4Subjects.mat');
load(filename);
if isempty(intersect(subjectNames,subjectNamesA))
    disp('details not saved for given subjects. Use saveRFdetails.m to save them first')
else
    [a,~,ind]= intersect(subjectNamesA,subjectNames,'stable');
%     disp(['Subjects available:',a])
    rfStatsDeg = rfStatsDeg(ind); 
    LFPElectrodeList =  LFPElectrodeList(ind);
    EcogElectrodeList = EcogElectrodeList(ind);
    subjectNames = subjectNames(ind);   
end
rfStatsPix=[];  % calculate separately

end