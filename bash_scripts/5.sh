 
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
  "bpd" "bpdz" "bpe" "bper" "bpet" "bpeu" "bpf" "bpg" "bph" "bpk" "bpl" "bpla" "bpm" "bpq" "elc" "eld" "elf" "elg" "elh" "eli" "ell" "eln" "elo" "elp" "elq" "elr" "elu" "elw" "ema" "eme" "emr" "hhs" "hht" "hhu" "hhy" "hhz" "hia" "hic" "hih" "hik" "hir" "hiw" "hix" "abv" "acal" "acan" "acav" "acc" "acd" "acep" "achb" "achr" "aci" "acia" "acid" "acii" "acin" "acio" "ack" "acn" "acom" "acp" "acr" "acra" "rop" "rpd" "rpe" "rpg" "rph" "rpha" "rpk" "rpl" "rpla" "rpn" "rpo" "rpp" "rpq" "rpr" "rps" "rpt" "rpv" "amar" "amb" "amed" "amf" "amg" "amil" "amim" "amk" "vex" "vff" "vfi" "vfl" "vfu" "vga" "vgo" "vha" "vhr" "vin" "vit" "vma" "vmi" "vmo" "vna" "vnl" "vok" "vow" "vpa" "vpb" "lmd" "lmf" "lmir" "lmn" "lmo" "lmob" "lmoc" "lmod" "lmoe" "lmog" "lmoj" "lmol" "lmom" "lmoq" 
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