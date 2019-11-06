 
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
  "kko" "kge" "ksd" "kpd" "mmw" "mme" "mpc" "tol" "tor" "oai" "mars" "bsan" "ncu" "nik" "gsn" "rfo" "ome" "aha" "ahy" "ahd" "ahr" "ahp" "ahj" "ahh" "ahi" "aaj" "asa" "aeo" "avr" "avo" "amed" "asr" "adh" "acav" "aem" "aea" "arv" "aes" "tau" "oce" "ocm" "opf" "zdf" "dno" "chj" "gap" "fpp" "sdf" "sok" "gbi" "slim" "sva" "acii" "saln" "tbn" "seds" "tsn" "thin" "tho" "pspi" "eof" "rma" "vok" "ebh" "bci" "bcib" "bcig" "gpb" "enm" "nme" "nmp" "nmh" "nmd" "nmm" "nms" "nmq" "nmz" "nma" "nmw" "nmx" "nmc" "nmn" "nmt" "nmi" "ngo" "ngk" "nla" "nel" "nwe" "nsi" "nmj" "nei" "nek" "nfv" "nsf" "salv" "kki" "vff" "vit" "ecor" "eex" "smur" "nba" "cvi" "cvc" "chro" "chri" "crz" "chrb" "chrm" "iod" "lhk" "pse" "jeu" "aql" "amah" "aqs" "rso" "rsc" "rsl" "rsn" "rsm" "rse" 
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