 
#!/bin/bash

#==============================================================
#**********
#**********  IGOR BRANDAO 2019 ©
#**********
#**********  PGM: TRIGGER SCRIPT TO DOWNLOAD ALL KEGG PATHWAYS
#**********
#**********  CLIENT: MASTER'S THESYS
#**********
#**********  VERSION:	1.0
#**********
#**********  APR/2019 - Creation
#**********
#==============================================================

#==============================================================
# GLOBAL VARIABLES
#==============================================================

NUMBER_OF_SCRIPTS=5
FOLDER_NAME='./bash_scripts'
SCRIPT_NAME='getKEGGData'

#==============================================================
# MAIN FLOW
#==============================================================

# Iterate over the organism list
for((i=2; i <= "${NUMBER_OF_SCRIPTS}"; i++)); do 
    bash ${FOLDER_NAME}'/'${SCRIPT_NAME}${i}'.sh' & 
done
