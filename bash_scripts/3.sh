 
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
  "tcl" "tcm" "tco" "tcx" "tdi" "tdl" "tdn" "tea" "tee" "teg" "ten" "oaq" "oar" "oat" "obb" "obo" "obr" "obt" "oca" "ocd" "ocg" "ock" "ocm" "ocn" "oco" "oct" "odi" "lap" "laq" "laqu" "lar" "las" "lat" "lau" "law" "lax" "lay" "laz" "lba" "lbc" "lbh" "lbi" "lbj" "lbk" "lbl" "lbo" "lbr" "lbu" "lby" "lca" "lcb" "lcc" "lcd" "lce" "pvu" "pvx" "pwo" "pxl" "pxy" "pya" "pyc" "pyg" "pyn" "pyo" "pys" "pzh" "qsu" "rab" "rac" "rae" "raf" "rak" "ram" "cmp" "cmw" "cna" "cne" "coc" "coh" "cohn" "col" "cola" "asoc" "asol" "asp" "asr" "ass" "asu" "asv" "asw" "asy" "asz" "atd" "atf" "ath" "ati" "atr" "ats" "fcm" "fco" "fcy" "fek" "fer" "ffa" "fgg" "fgl" "fgo" "fgu" "fhw" "fil" "fin" "fiy" "fki" "fla" "fli" "flm" "fln" "fmi" "fmo" "hpd" "hpf" "hph" "hpit" "hpn" "hpr" 
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