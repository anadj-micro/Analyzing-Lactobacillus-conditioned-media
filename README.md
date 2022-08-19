# Analyzing-Lactobacillus-conditioned-media
Data and scripts provided here are meant for comparing the secretomes of L. murinus, L. rhamnosus, and E. coli (used as a negative control). To do this we grew these strains in BHI media, collected their supernatants, extracted metabolites, derivatized them and analyzed on Agilent 8890 GC coupled to Agilent 5977B mass selective detector. Metabolomics data have been deposited to the EMBL-EBI  MetaboLights database (DOI: 10.1093/nar/gkz1019, PMID:31691833) with the  identifier MTBLS5733. The complete dataset can be accessed here https://www.ebi.ac.uk/metabolights/MTBLS5733. Obtained peak areas, corresponding to different metabolites, were subjected to the analysis described here. 
Download the raw data from MetaboLights depository (link above). Download the chromatography-master directory from here. The raw data should be in the folder outside the folder where script and chromatography-master directory are located. This directory should be titled "MRE_metabolomics". tblSamplesSubmission.txt contains metadata necessary for processing the peaks. Download it to the same directory from which you are running the script (scAnalyzeLactobacillusSecretome.m)
