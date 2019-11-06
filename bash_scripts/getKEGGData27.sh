 
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
  "sje" "apak" "aqt" "psyh" "psyo" "prd" "stea" "bse" "sau" "sav" "saw" "sah" "saj" "sam" "sas" "sar" "sac" "sax" "saa" "sao" "sae" "sad" "suu" "suv" "sue" "suj" "suk" "suc" "sut" "suq" "suz" "sud" "sux" "suw" "sug" "suf" "saua" "saue" "saun" "saus" "sauu" "saug" "sauz" "saut" "sauj" "sauk" "sauq" "sauv" "sauw" "saux" "sauy" "sauf" "sab" "suy" "saub" "saum" "sauc" "saur" "saui" "saud" "sams" "suh" "sep" "ser" "sepp" "seps" "sha" "shh" "ssp" "sca" "slg" "sln" "ssd" "sdt" "swa" "spas" "sxy" "sxl" "sxo" "shu" "scap" "ssch" "sscz" "sagq" "seqo" "ssif" "scv" "spet" "slz" "scoh" "sscu" "snl" "skl" "sfq" "shom" "smus" "mcl" "mcak" "macr" "shv" "sbac" "lmo" "lmn" "lmy" "lmt" "lmoc" "lmoe" "lmob" "lmod" "lmow" "lmoq" "lmr" "lmom" "lmf" "lmc" "lmog" "lmp" "lmol" "lmoj" "lmoz" 
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