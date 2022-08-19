%% Ana Djukovic, Ph.D., Memorial Sloan-Kettering Cancer Center, Joao Xavier lab
% This script is intended to compare secretomes of L. murinus, L. rhamnosus,
% and E. coli (used as a negative control). To do this we grew these strains in BHI
% media, collected their supernatants and analyzed on Agilent 8890 GC coupled to
% Agilent 5977B mass selective detector. Peak areas were subjected to
% following analysis. Metabolomics data have been deposited to the EMBL-EBI 
% MetaboLights database (DOI: 10.1093/nar/gkz1019, PMID:31691833) with the 
% identifier MTBLS5733. The complete dataset can be accessed here
% https://www.ebi.ac.uk/metabolights/MTBLS5733. In this script the raw data
% is located in the one-level-up directory named "MRE_metabolomics"-
% change accordingly to your directory organization. chromatography-master
% directory is necessary for running the script, download it as well.

% Last updated: Aug. 19, 2022.

% load samples table
tblSamples = readtable('tblSamplesSubmission.txt', 'Delimiter', '\t');
tblSamples.sampleId = categorical(tblSamples.sampleId);

%% add columns that define replicates and groups
for i = 1:height(tblSamples)
    s = split(tblSamples.sample{i}, '_');
    s = double(uint8(s{1}(end)));
    s = char(s);
    tblSamples.replicate{i} = s;
end
tblSamples.groups = categorical(join(tblSamples{:, {'sampletype' 'cell_and_color'}}));


%% We use the raw data (the GCMS chromatograms) to find compounds
% in media conditioned by Lactobacillus murinus or L. rhamnosus compared 
% to the fresh media. Include E.coli as a control.

% load GCMS raw data
dataDir = '../MRE_metabolomics'; %change to the directory where you keep the raw data

% add the chromatography master package
addpath(genpath('chromatography-master'));

% import GCMS files
h = waitbar(0,'Loading files...');
for i = 1:height(tblSamples)
    sId = strsplit(tblSamples.sample{i}, '_');
    sId = sId{1};
    fn = sprintf('%s/%s_split.D', dataDir, sId);
    data.sample{i} = importAgilentAndResample(fn, sId);
    waitbar(i/height(tblSamples),h)
end
close(h);

%% build the data table that contains information about peak areas for each m/z ("mz" column) and retention time ("time" column) 
mz = data.sample{1}.mz;
time = data.sample{1}.time;

[mzMesh, timeMesh] = meshgrid(mz, time);
tlbSpectra = table();
tlbSpectra.mz = mzMesh(:);
tlbSpectra.time = timeMesh(:);

for i = 1:height(tblSamples)
    tlbSpectra{:, string(tblSamples.sampleId(i))} = data.sample{i}.xic(:);
end

%% stack the data
tlbSpectraStacked = stack(tlbSpectra, string(tblSamples.sampleId),...
    'NewDataVariableName','xic',...
    'IndexVariableName','sample');

%% make the Total Ion Chromatogram (TIC) for each and all the samples
tlbTicOfSamples = tlbSpectraStacked;
tlbTicOfSamples.mz = [];
tlbTicOfSamples = unstack(tlbTicOfSamples,'xic', 'sample');

tlbTicOfAllSamples = tlbTicOfSamples(:, "time");
tlbTicOfAllSamples.tic = sum(tlbTicOfSamples{:, 2:end}, 2);

%% find the peaks in the TIC
[Peaklist, PFWHH, PExt] = mspeaks(tlbTicOfAllSamples.time, tlbTicOfAllSamples.tic);

%% integrate TIC for each sample and  each peak
tblPeaks = table(Peaklist(:,1), Peaklist(:,2), PExt(:, 1), PExt(:, 2));
tblPeaks.Properties.VariableNames = {'rt', 'ticSum', 'rtStart', 'rtEnd'};
for i = 1:height(tblPeaks)
    idx = tlbTicOfSamples.time >= tblPeaks.rtStart(i) &...
        tlbTicOfSamples.time < tblPeaks.rtEnd(i);
    % integrate peak
    tblPeaks{i, string(tblSamples.sampleId)} = sum(tlbTicOfSamples{idx, string(tblSamples.sampleId)}, 1);
end

%% stack the tblPeaks
tblPeaksStacked = stack(tblPeaks, string(tblSamples.sampleId),...
    "NewDataVariableName", 'peakArea',...
    "IndexVariableName", 'sampleId');

tblPeaksStacked = innerjoin(tblPeaksStacked,...
    tblSamples(:, {'sampleId' 'sampletype' 'cell_and_color' 'groups' 'replicate'}));


%% do anova for each peak to determine if it changes among different supernatants/un-inoculated media
h = waitbar(0,'Anova for each peak...');
for i = 1:height(tblPeaks)
    idx = tblPeaksStacked.rt == tblPeaks.rt(i);
    % anova to determine if metabolite is relevant
    mdl = fitlm(tblPeaksStacked(idx, :), 'peakArea ~ groups');
    tblAnova = anova(mdl);
    tblPeaks.anovaPValue(i) = tblAnova.pValue(1);
    % calculate the log2FC in samples using BHI (un-inoculated media) as a reference
    tblT0 = tblPeaksStacked(idx, :);
    tblT0(~contains(tblT0.sampletype, 'T0'), :) = [];
    cats = {'BHI',...
        'E_coli', 'L_murinus', 'L_rhamnosus'};
    tblT0.cell_and_color = categorical(tblT0.cell_and_color, {'BHI',...
        'E_coli', 'L_murinus', 'L_rhamnosus'});
    tblT0.peakArea = log2(tblT0.peakArea);
    mdl2 = fitlme(tblT0, 'peakArea ~ cell_and_color + (1|replicate)');
    for j = 2:length(cats)
        tblPeaks{i, [cats{j} '_Log2Fc']} = mdl2.Coefficients.Estimate(j);
        tblPeaks{i, [cats{j} '_pValue']} = mdl2.Coefficients.pValue(j);
    end
    waitbar(i/height(tblPeaks),h)
end
close(h);

writetable(tblPeaks, 'tblPeaks.txt', 'Delimiter', '\t')

%% scatter the FC to find metabolites that changed
%pCutOff = 0.05 / height(tblPeaks);
pCutOff = 0.05 ;

% pFlag has value 1 if peak changed in rhamnosus or murinus, and value 2
% if it changed in both
pFlag = (tblPeaks.L_murinus_pValue < pCutOff) +...
    (tblPeaks.L_rhamnosus_pValue < pCutOff);

figure(1)
subplot(1, 2, 1)
scatter(tblPeaks.L_murinus_Log2Fc, tblPeaks.L_rhamnosus_Log2Fc,...
    [], pFlag, 'filled')
xlabel('log_2(FC) L. murinus');
ylabel('log_2(FC) L. rhamnosus');
grid on
title('Metabolite peaks changed in Lactobacillae conditioned media')

subplot(1, 2, 2)
scatter(tblPeaks.L_murinus_Log2Fc, tblPeaks.E_coli_Log2Fc,...
    [], pFlag, 'filled')
xlabel('log_2(FC) L. murinus');
ylabel('log_2(FC) E. coli');
grid on
title('Compared with E. coli supernatant')

%% bar plot of fold changes

figure(2)
barh([tblPeaks.L_murinus_Log2Fc(pFlag>0),...
    tblPeaks.L_rhamnosus_Log2Fc(pFlag>0),...
    tblPeaks.E_coli_Log2Fc(pFlag>0)])
set(gca, 'YTickLabel', tblPeaks.rt(pFlag>0))
xlabel('log_2(FC) [from fitlm effect sizes]')
ylabel('Peak retention time [min]')
legend('L.murinus', 'L.rhamnosus', 'E.coli')
title('Peaks significantly changed')

