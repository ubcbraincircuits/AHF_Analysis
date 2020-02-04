# AHF
MATLAB code for 2020 AHF paper

MAIN scripts are the primary processing pipelines. They associate behavior data, text files, and brain data; have various automated and manual quality control steps; perform preprocessing and time alignment; and save processed outputs as .mat files. These scripts call the FARM functions. Each MAIN script targets a specific group of mice and/or range of days.

STEP2 script is the analysis pipeline for figures 9 and 10 in the manuscript. STEP2 calls dffPxlDecoder, shadedErrorBar, and violinplot. Additional helper functions are listed at the bottom of the script. STEP2 uses the outputs from the MAIN scripts.
