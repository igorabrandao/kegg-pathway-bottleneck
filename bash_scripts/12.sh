 
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
  "nao" "nap" "nat" "nau" "nba" "nbg" "nmj" "nmm" "nmn" "nmp" "nmq" "nms" "nmt" "nmu" "nmw" "nmx" "nob" "noc" "noj" "nok" "nom" "nou" "npa" "npe" "npl" "otm" "ots" "ott" "our" "ovi" "pab" "pac" "paca" "pace" "pach" "pacn" "pacr" "parc" "paro" "parr" "pars" "part" "paru" "pary" "pat" "pau" "paur" "pavl" "pay" "pazo" "pbar" "pbb" "pbc" "pbh" "pfa" "pfae" "pfb" "pfc" "pfd" "pfe" "pff" "pfg" "pfh" "pfi" "pfj" "pfk" "pfl" "pfn" "pfo" "pfp" "pfr" "pfre" "pfs" "pfu" "pfv" "pll" "pln" "plo" "plq" "plu" "plum" "plv" "plw" "plx" "ply" "pmai" "pman" "pmar" "pmib" "psv" "psw" "pswu" "psx" "psy" "psyg" "psyr" "psz" "pta" "ptc" "pte" "ptep" "pth" "pti" "ptl" "sedi" "seds" "see" "seeh" "seen" "sef" "seg" "sega" "seh" "sei" "sek" "sel" "sele" "selo" "selt" "sena" "send" 
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