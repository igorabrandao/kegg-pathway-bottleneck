 
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
  "pand" "plg" "hyf" "lmir" "mcys" "bpe" "bpc" "bper" "bpet" "bpeu" "bpa" "bpar" "bbr" "bbm" "bbh" "bbx" "bpt" "bav" "bho" "bhm" "bhz" "btrm" "bbro" "bfz" "bpdz" "boh" "bgm" "boz" "boj" "axy" "axo" "axn" "axx" "adt" "ais" "asw" "achr" "achb" "teq" "tea" "teg" "tas" "tat" "put" "pus" "aka" "amim" "cdn" "bpsi" "afa" "afq" "aaqu" "phn" "odi" "our" "pig" "pacr" "kgy" "rfr" "rsb" "rac" "rhy" "rhf" "rhg" "pol" "pna" "pos" "aav" "ajs" "dia" "aaa" "ack" "acra" "acid" "acip" "acin" "acis" "acio" "vei" "dac" "del" "dts" "dhk" "vap" "vpe" "vpd" "vaa" "vbo" "vam" "ctt" "ctes" "cke" "cser" "cof" "adn" "adk" "rta"
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