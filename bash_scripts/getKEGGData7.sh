 
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
  "ddd" "dze" "ddc" "dzc" "dso" "ced" "dfn" "ddq" "daq" "dic" "bgj" "brb" "bng" "lbq" "eta" "epy" "epr" "eam" "eay" "ebi" "erj" "ege" "epe" "ehd" "buc" "bap" "bau" "baw" "bajc" "bua" "bup" "bak" "buh" "bapf" "bapg" "bapu" "bapw" "bas" "bab" "bcc" "baj" "baph" "wbr" "wgl" "pam" "plf" "paj" "paq" "pva" "pagg" "pao" "kln" "pant" "panp" "hhs" "pck" "pagc" "pstw" "palh" "pans" "pgz" "pcd" "tci" "tpty" "plu" "plum" "pay" "ptt" "pmr" "pmib" "pvl" "pvg" "phau" "xbo" "xbv" "xne" "xnm" "xdo" "xpo" "xho" "psi" "psx" "psta" "prg" "pala" "phei" "prq" "prj" "mmk" "asy" "aen" "ans" "eic"
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