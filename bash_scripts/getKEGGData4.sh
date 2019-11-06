 
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
  "ecoo" "ecoh" "ecg" "eok" "elr" "eso" "esm" "esl" "ecw" "elh" "eun" "ecc" "ecp" "eci" "ecv" "ecx" "ecm" "ecy" "ecr" "ecq" "eck" "ect" "eoc" "eum" "ecz" "elo" "eln" "ese" "ecl" "ebr" "ebd" "eko" "ekf" "eab" "edh" "edj" "eih" "ena" "elu" "elw" "ell" "elc" "eld" "elp" "ebl" "ebe" "elf" "ecoa" "ecol" "ecoi" "ecoj" "ecos" "efe" "eal" "ema" "sty" "stt" "sex" "sent" "stm" "seo" "sev" "sey" "sem" "sej" "seb" "sef" "setu" "setc" "senr" "send" "seni" "seen" "spt" "sek" "spq" "sei" "sec" "seh" "shb" "senh" "seeh" "see" "senn" "sew" "sea" "sens" "sed" "seg" "sel" "sega" "set" "sena" "seno" "senv" "senq" "senl"
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