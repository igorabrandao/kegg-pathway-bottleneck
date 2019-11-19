 
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
  "aad" "aae" "aaf" "aah" "aaj" "aal" "aalg" "aalt" "aamy" "aan" "aao" "aap" "aaqu" "aat" "aaus" "aaw" "aay" "mcg" "mch" "mche" "mcj" "mcl" "mcn" "mcs" "mct" "mcw" "mcys" "eal" "eam" "ean" "ear" "eas" "eat" "eau" "eay" "eba" "ebc" "ebd" "ebe" "ebf" "ebh" "ebi" "ebl" "ebr" "ebs" "ebu" "phl" "pho" "php" "phq" "phr" "phu" "phz" "pia" "pib" "pic" "pif" "pig" "pih" "pin" "pjd" "pkb" "pkc" "sbs" "sbt" "sbw" "sca" "scap" "sch" "scj" "sck" "scl" "scla" "sclo" "scm" "scoh" "scon" "scos" "scou" "ano" "app" "att" "axe" "babr" "bah" "bapw" "bbo" "bcen" "bgf" "bho" "blb" "ccav" "cci" "ccm" "cco" "ccoc" "ccon" "gho" "gjf" "gka" "gli" "glo" "gmc" "gme" "gmo" "gmx" "gni" "goh" "gor" "got" "gox" "goy" "gpb" "gpi" "tbw" "tca" "tcb" "tcc" "tci" 
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