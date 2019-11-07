 
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
  "bml" "bmn" "bmal" "bmae" "bmaq" "bmai" "bmaf" "bmaz" "bmab" "bps" "bpm" "bpl" "bpd" "bpr" "bpse" "bpsm" "bpsu" "bpsd" "bpz" "bpq" "bpk" "bpsh" "bpsa" "bpso" "but" "bte" "btq" "btj" "btz" "btd" "btv" "bthe" "bthm" "btha" "bthl" "bok" "boc" "buu" "bvi" "bve" "bur" "bcn" "bch" "bcm" "bcj" "bcen" "bcew" "bceo" "bam" "bac" "bmu" "bmj" "bmk" "bmul" "bct" "bced" "bcep" "bdl" "bpyr" "bcon" "bub" "bdf" "blat" "btei" "bsem" "bpsl" "bmec" "bstg" "bstl" "bgl" "bgu" "bug" "bgf" "bgd" "bgo" "byi" "buk" "buo" "bue" "bul" "buq" "bgp" "bpla" "bud" "buz" "bum" "bui" "bxe" "bxb" "bph" "bge" "bpx" "bpy" "bfn" "bcai" "pspw" "para"
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