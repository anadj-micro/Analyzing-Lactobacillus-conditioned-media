function sample=importAgilentAndResample(fn, sampleId)

%importAgilentAndResample: imports Agilent data and resamples it to a fixed
% grid or retention time to mz.

% import raw data using the ImportAgilent function from
% chromatography-master 
hiRes = ImportAgilent(fn, 'precision', 3);

% copy but downsample to a common grid of rt and mz
sample      = hiRes;
sample.time = (6:0.01:30)';
sample.mz   = (50:599); 

% subsample each eic by adding all the abundances detected
sample.xic = zeros(length(hiRes.time), length(sample.mz));
for i = 1:(length(sample.mz)-1)
    j = (hiRes.mz >= sample.mz(i)) & (hiRes.mz < sample.mz(i+1));
    sample.xic(:, i) = sum(hiRes.xic(:, j), 2);
end
% last bin
j = (hiRes.mz >= sample.mz(end));
sample.xic(:, end) = sum(hiRes.xic(:, j), 2);

% subsample the residence times using linear interpolation
% maybe another way would be subsample with loess
sample.xic = interp1(hiRes.time, sample.xic, sample.time);

% % recalculate the tic
% sample.tic = sum(sample.xic, 2);
% 
% % add the other atttributes
% sample.sampleId = sampleId;
% sample.sampleType = [cellType '_' labeling];
% sample.cellType = cellType;
% sample.labeling = labeling;







