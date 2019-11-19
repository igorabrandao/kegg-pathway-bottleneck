 
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
  "bsem" "bsf" "bsg" "bsh" "bsi" "bsj" "bsk" "bsl" "bsm" "bsn" "bso" "bsp" "bsq" "bsr" "bss" "bst" "bstl" "bsto" "bsu" "bsuc" "bsui" "bsul" "bsup" "bsus" "bsut" "bsuv" "bsv" "bsw" "bsx" "bsy" "bsz" "btab" "beo" "beq" "bex" "bfm" "bfx" "bfz" "bgd" "bge" "ecm" "ecoa" "ecog" "ecoh" "ecoi" "ecoj" "ecol" "ecor" "ecos" "ecp" "ecq" "ecr" "ect" "ecv" "ecw" "ecx" "ecy" "cfx" "cfz" "cgi" "cgr" "cgrn" "cgw" "cha" "sod" "soe" "soi" "sok" "son" "soo" "sot" "sox" "soz" "spa" "span" "spas" "spat" "spb" "spc" "pcr" "pcre" "pcs" "pcx" "pcy" "pcz" "pda" "pdag" "pdam" "pdg" "pdh" "pdi" "pdp" "pdr" "pds" "pdu" "pea" "pef" "pel" "pen" "pep" "peq" "pes" "pet" "hpya" "hpyb" "hpyc" "hpyd" "hpye" "hpyf" "hpyg" "hpyh" "hpyi" "hpyj" "hpyk" "hpym" "hpyo" "hpyr" "hpys" "hpyu" "hpz" 
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