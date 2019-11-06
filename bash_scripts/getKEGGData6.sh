 
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
  "kie" "icp" "lax" "lei" "leh" "lee" "laz" "lef" "lni" "lew" "lpv" "buf" "sbw" "den" "hed" "ged" "cmik" "ppet" "mety" "ahn" "izh" "ebf" "ebc" "ebu" "psts" "ype" "ypk" "yph" "ypa" "ypn" "ypm" "ypp" "ypg" "ypz" "ypt" "ypd" "ypx" "ypw" "ypj" "ypv" "ypl" "yps" "ypo" "ypi" "ypy" "ypb" "ypq" "ypu" "ypr" "ypc" "ypf" "yen" "yep" "yey" "yel" "yew" "yet" "yef" "yee" "ysi" "yal" "yfr" "yin" "ykr" "yro" "yru" "yrb" "yak" "yma" "yhi" "spe" "srr" "srl" "sry" "sply" "srs" "sra" "ssz" "smaf" "smw" "smar" "smac" "slq" "serf" "sers" "sfw" "sfg"
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