# AHF
MATLAB code for 2020 AHF paper

MAIN scripts are the primary processing pipelines. They associate behavior data, text files, and brain data; have various automated and manual quality control steps; perform preprocessing and time alignment; and save processed outputs as .mat files. These scripts call the FARM functions. Each MAIN script targets a specific group of mice and/or range of days.

STEP2 script is the analysis pipeline for figures 9 and 10 in the manuscript. STEP2 calls dffPxlDecoder, shadedErrorBar, and violinplot. Additional helper functions are listed at the bottom of the script. STEP2 uses the outputs from the MAIN scripts.

#AutoHeadFix 2019 Female Cage Readme 
All data was stored in MariaDB(v10.1.38-MariaDB-0+deb9u1 Raspbian 9.0) and accessed via SQL queries
Data was then analyzed in python (v3.7.4) including modules matplotlib(v3.1.1), numpy(1.17.0), pandas(v0.25.4), 
pymysql(v0.9.3), and pytz(v2019.2).
Intruction to Access/Run IPython script
1. Import the database dump file (AHF2019_Females.sql, available at the paper's subsequent scholarverse page) to a suitable computer with MySQL database.(Works with MySQL 8.0)
2. Follow the instructions in cell #2 of IPython notebook (AHF_2019_fucntion_analysis.ipynb) and fill out the appropiate credientials (which is based on your own db login user and password) to run the script.
Note: The IPython notebook will generate and saves subsequent data in csv or svg formate in the same folder
