 
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
  "clr" "clm" "clq" "cln" "cll" "ccol" "ccc" "ccq" "ccf" "ccy" "ccoi" "ccof" "ccoo" "caj" "cis" "cvo" "cpel" "camr" "csm" "csf" "cgra" "cure" "chyo" "chv" "cspf" "cpin" "ccun" "clx" "cavi" "chw" "camz" "camy" "abu" "abt" "abl" "ant" "ask" "ahs" "amyt" "amar" "arc" "alp" "sdl" "sba" "smul" "shal" "suls" "sulj" "hyo" "nsa" "nis" "sun" "slh" "nam" "nap" "cmed" "cpaf" "gsu" "gsk" "gme" "gur" "glo" "gbm" "geo" "gem" "geb" "gpi" "gao" "gsb" "pca" "ppd" "pace" "pef" "des" "deu" "dvu" "dvl" "dvm" "dvg" "dde" "dds" "ddn" "dma" "dsa" "dhy" "dgg" "dfi" "dpg" "def" "dtr" "dfl" "daf" "das" "dpi" "dej" "pprf" "lip" "lir" "dba" "doa" "drt" "dps" "dak" "dpr" "deo" "dsf" "dol" "dml" "dal" "dat" "dto" "ade" "acp" "afw" "ank" "mxa" "mfu" "msd" "mym" "mmas" "ccx"
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