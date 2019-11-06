 
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
  "acet" "abg" "kba" "rgi" "ros" "rmuc" "nch" "coq" "shum" "ntn" "neh" "ssam" "swf" "rru" "rrf" "rce" "rpm" "mag" "mgy" "mgry" "magx" "magn" "azl" "ali" "abs" "abq" "abf" "ati" "ahu" "azt" "azm" "azz" "tmo" "thal" "efk" "txi" "thac" "tii" "magq" "hjo" "nao" "ncb" "fer" "pbr" "mgm" "pub" "pel" "apc" "apm" "ecog" "mai" "man" "pgv" "phr" "pstg" "apb" "afe" "afr" "acu" "acz" "afi" "afj" "maes" "mfn" "htl" "bsu" "bsr" "bsl" "bsh" "bsy" "bsut" "bsul" "bsus" "bss" "bst" "bso" "bsn" "bsq" "bsx" "bsp" "bli" "bld" "blh" "bay" "baq" "bya" "bamp" "baml" "bama" "bamn" "bamb" "bamt" "bamy" "bmp" "bao" "baz" "bql" "bxh" "bqy" "bami" "bamc" "bamf" "bae" "bvm" "bha" "ban" "bar" "bat" "bah" "bai" "bax" "bant" "banr" "bans" "banh" "banv" "bce" "bca" "bcz" "bcr" "bcb" "bcu" "bcg" 
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