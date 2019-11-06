 
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
  "psya" "psyy" "psyp" "acb" "abm" "aby" "abc" "abn" "abb" "abx" "abz" "abr" "abd" "abh" "abad" "abj" "abab" "abaj" "abaz" "abk" "abau" "abaa" "abw" "abal" "acc" "ano" "alc" "acal" "acd" "aci" "att" "aei" "ajo" "acw" "acv" "ahl" "ajn" "asol" "ala" "asj" "aid" "adv" "arj" "awu" "acum" "mct" "mcs" "mcat" "moi" "mos" "mbl" "mboi" "mbah" "son" "sdn" "sfr" "saz" "sbl" "sbm" "sbn" "sbp" "sbt" "sbs" "sbb" "slo" "spc" "shp" "sse" "spl" "she" "shm" "shn" "shw" "shl" "swd" "swp" "svo" "shf" "sja" "spsw" "sbj" "smav" "shew" "salg" "slj" "ilo" "ili" "ipi" "idi" "idt" "cps" "com" "coz" "colw" "cola" "cber" "cov" "lsd" "tht" "thap" "pha" "pat" "psm" "pseo" "pia" "pphe" "pbw" "prr" "ptn" "plz" "paln" "ppis" "pea" "pspo" "part" "ptu" "png" "ptd" "psen" "pdj" 
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