 
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
  "amg" "amk" "alt" "aal" "aaus" "asp" "asq" "aaw" "alr" "ale" "alz" "gag" "gni" "gps" "lal" "cate" "salh" "salm" "salk" "hmi" "pin" "psy" "fbl" "mvs" "mya" "mmaa" "cja" "ceb" "cell" "cek" "sde" "ttu" "saga" "spoi" "zal" "osg" "mthd" "micc" "maga" "mii" "hja" "cbu" "cbs" "cbd" "cbg" "cbc" "cea" "cey" "rvi" "alg" "lpn" "lph" "lpo" "lpu" "lpm" "lpf" "lpp" "lpc" "lpa" "lpe" "llo" "lfa" "lha" "lok" "lcd" "les" "lsh" "llg" "lib" "lgt" "tmc" "mca" "mmt" "mdn" "mdh" "mko" "metl" "mah" "mbur" "mpsy" "mmai" "ftu" "ftq" "ftf" "ftw" "ftr" "ftt" "ftg" "ftl" "fth" "fta" "fts" "fti" "fto" "ftc" "ftv" 
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