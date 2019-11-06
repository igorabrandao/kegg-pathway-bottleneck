 
#!/bin/bash

#==============================================================
#**********
#**********  IGOR BRANDAO 2019 ©
#**********
#**********  PGM: BASH SCRIPT TO DOWNLOAD ALL KEGG PATHWAYS
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

# Data type: ec/ko/MAP
DATA_TYPE="ko"

# Current directory
PATHPWD=$(pwd)

# Organisms list
declare -a ORG_LIST=(
  "cpae" "csb" "cah" "clt" "cbv" "csq" "cld" "cace" "cck" "cbut" "ctyk" "ceu" "ctae" "cfm" "cchv" "carg" "cdrk" "cia" "csep" "cdy" "amt" "aoe" "asf" "asm" "aso" "asb" "gfe" "hhw" "cale" "crs" "clo" "fsa" "cth" "ctx" "ccl" "hsc" "ruk" "cce" "css" "csd" "cthd" "esr" "esu" "ccel" "rbp" "fpla" "eha" "ral" "rch" "rum" "rus" "ruj" "fpr" "fpa" "capr" "bpb" "bfi" "bhu" "cle" "rho" "rix" "rim" "coo" "cct" "rob" "byl" "rto" "bhan" "blau" "bpro" "cpy" "lacy" "csh" "cso" "cbol" "bprl" "arf" "hsd" "cpro" "lua" "ehl" "pxv" "ere" "ert" "era" "lbw" "cdf" "pdc" "cdc" "cdl" "pdf" "eac" "cst" "faa" "psor" "pbq" "sth" "swo" "slp" "dsy" "dhd" "ddh" "ddl" "dmt" "drm" "dca" "dru" "dfg" "dae" "dku" "dgi" "pth" "dau" "tjr" "sgy" "dor" "dai" "dmi" "ded" "dec" "drs" "hmo" "eel" "elm" "emt" 
)

#==============================================================
# MAIN FLOW
#==============================================================

echo "<<< INITIALIZING KEGG KGML FILES DOWNLOAD..."

# Exit the bash_scripts folder
cd '..'

# Iterate over the organism list
for i in "${!ORG_LIST[@]}"; do
	echo "***************"
	echo ${ORG_LIST[$i]}
	echo "***************"

	# Generate current organism folder name
	org_folder="./output/kgml/${ORG_LIST[$i]}"

  # Check if the browser folder exists
  if [ ! -d $org_folder ]; then
      mkdir -p $org_folder;
  fi

  # Enter into the org folder data
  cd $org_folder

	# Perform the KGML file download
	curl -s http://rest.kegg.jp/list/pathway/${ORG_LIST[$i]} | awk '{split($1,a,":"); print "curl -s http://rest.kegg.jp/get/"a[2]"/kgml -o ./"a[2]".xml"}' | bash

  # Exit the org folder data
  cd '..'  
  cd '..'
  cd '..'

done

echo "<<< PROCESS FINISHED SUCCESSFULLY!"

#==============================================================