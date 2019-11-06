 
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
  "senj" "seec" "seeb" "seep" "senb" "sene" "senc" "ses" "sbg" "sbz" "sbv" "salz" "sfl" "sfx" "sfv" "sfe" "sfn" "sfs" "sft" "ssn" "sbo" "sbc" "sdy" "sdz" "shq" "ent" "enc" "eno" "eclo" "enl" "eclg" "ecle" "ecln" "ecli" "eclx" "ecly" "eclz" "ehm" "exf" "ecla" "eclc" "eau" "ekb" "eec" "elg" "ecan" "ern" "ecls" "echg" "eas" "enr" "enx" "enf" "ebg" "esa" "csk" "csz" "csj" "ccon" "cdm" "csi" "cmj" "cui" "cmw" "ctu" "kpn" "kpu" "kpm" "kpp" "kph" "kpz" "kpv" "kpw" "kpy" "kpg" "kpc" "kpq" "kpt" "kpe" "kpo" "kpr" "kpj" "kpi" "kpa" "kps" "kpx" "kpb" "kpne" "kpnu" "kpnk" "kva" "kpk" "kvd" "kvq" "kox" "koe" "koy" "kom" "kmi" "kok" "koc" "kqu" "eae" "ear" "kqv" "kll" "klw" "cko" "cro" "cfd" "cbra" "cwe" "cyo" "cpot" "cama" "caf" "cif" "cfar" "cir" "cie" "cpar"
)

#==============================================================
# MAIN FLOW
#==============================================================

echo "<<< INITIALIZING KEGG KGML FILES DOWNLOAD..."

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