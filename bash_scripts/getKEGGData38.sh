 
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
  "ltr" "agg" "art" "arr" "arm" "arl" "are" "aaq" "arw" "arh" "ary" "arz" "aru" "arq" "arn" "arx" "acry" "arth" "artp" "ari" "aau" "ach" "apn" "psul" "psni" "aai" "gar" "gcr" "glu" "rsa" "krh" "kpl" "kfv" "kii" "krs" "mlu" "mick" "rmu" "rdn" "raj" "satk" "nae" "aul" "cig" "bcv" "bfa" "brx" "brv" "bgg" "brz" "bsau" "brr" "dva" "djj" "jde" "kse" "dni" "day" "lmoi" "xce" "iva" "ido" "cet" "cceu" "xya" "xyl" "ske" "cfl" "cfi" "cga" "cez" "celz" "oek" "ica" "ars" "serj" "serw" "jte" "jli" "orn" "orz" "bly" "blin" "bri" "dco" "gez" "pac" "pak" "pav" "pax" "paz" "paw" "pad" "pcn" "pacc" "pach" "pacn" "cacn" "pra" "cgrn" "pfr" "pfre" "prl" "pacd" "ppc" "pbo" "aaci" "acij" "aji" "mph" "mik" "tfl" "tfa" "tes" "tez" "nca" "ndk" "noy" "noi" "kfl" "psim" "aer" "aez" "aeb" 
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