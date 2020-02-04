AutoHeadFix 2019 Female Cage Readme 
All data was stored in MariaDB(v10.1.38-MariaDB-0+deb9u1 Raspbian 9.0) and accessed via SQL queries
Data was then analyzed in python (v3.7.4) including modules matplotlib(v3.1.1), numpy(1.17.0), pandas(v0.25.4), 
pymysql(v0.9.3), and pytz(v2019.2).
##############################################
##############################################
##############################################
Intruction to Access/Run IPython script
1. Import the database dump file (AHF2019_Females.sql) to a suitable computer with MySQL database.(Works with MySQL 8.0)
2. Follow the instructions in cell #2 of IPython notebook (AHF_2019_fucntion_analysis.ipynb) and fill out the appropiate credientials to run the script.
Note: The IPython notebook will generate and saves subsequent data in csv or svg formate in the same folder