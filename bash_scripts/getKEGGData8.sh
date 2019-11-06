 
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
  "hiw" "hic" "hix" "hpr" "hdu" "hay" "hpit" "hhz" "hap" "hpaz" "hpas" "hpak" "gle" "hso" "hsm" "pmu" "pmv" "pul" "pmp" "pmul" "pdag" "paet" "msu" "mht" "mhq" "mhat" "mhx" "mhae" "mham" "mhao" "mhal" "mhaq" "mhay" "mvr" "mvi" "mvg" "mve" "apl" "apj" "apa" "asu" "asi" "ass" "aeu" "apor" "aio" "aap" "aaz" "aat" "aao" "aan" "aah" "aacn" "aact" "aseg" "gan" "bto" "btre" "btrh" "btra" "apag" "rpne" "xfa" "xft" "xfm" "xfn" "xff" "xfl" "xfs" "xfh" "xtw" "xcc" "xcb" "xca" "xcp" "xcv" "xax" "xac" "xci" "xct" "xcj" "xcu" "xcn" "xcw" "xcr" "xcm" "xcf" "xfu" "xao" "xoo" "xom" "xop" "xoy" "xor" "xoz" 
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