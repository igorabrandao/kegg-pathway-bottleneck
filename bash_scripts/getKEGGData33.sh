 
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
  "awo" "ova" "obj" "tmr" "thef" "say" "sap" "sthr" "cthm" "cmiu" "ibu" "mdv" "amij" "euu" "bprm" "bprs" "cbar" "tte" "tex" "thx" "tpd" "tit" "tmt" "tbo" "twi" "tki" "chy" "tep" "tae" "mta" "adg" "tpz" "csc" "ate" "cob" "chd" "cow" "cki" "ckn" "clc" "ccha" "toc" "ttm" "tto" "txy" "tsh" "tnr" "taci" "mas" "nth" "hor" "has" "hpk" "hals" "aar" "hhl" "aft" "fma" "apr" "pmic" "ped" "phar" "cad" "spoa" "vpr" "vat" "vrm" "vdn" "med" "mhw" "meg" "dpn" "ssg" "sri" "sele" "selo" "selt" "mhg" "puf" "pft" "mana" "sted" "afn" "ain" "pfac" "erh" "ers" "erl" "eri" "euc" "tur" "fro" "erb" "ebm" "lpil" "cpo" "mge" "mgu" "mgc" "mgq" "mgx" "mpn" "mpm" "mpj" "mpb" "mpu" "mpe" "mga" "mgh" "mgf" "mgn" "mgs" "mgt" "mgv" "mgw" "mgac" "mgan" "mgnc" "mgz" "mmy" "mmym" "mmyi" "mml" "mcp" 
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