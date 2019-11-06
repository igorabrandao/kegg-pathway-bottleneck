 
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
  "bcq" "bcx" "bal" "bnc" "bcf" "bcer" "bcef" "bcy" "btk" "btl" "btb" "btt" "bthr" "bthi" "btc" "btf" "btm" "btg" "bti" "btn" "btht" "bthu" "btw" "bthy" "bwe" "bww" "bmyo" "bty" "bmyc" "bby" "bwd" "bcl" "bpu" "bpum" "bpus" "bpf" "bmq" "bmd" "bmh" "bmeg" "bco" "bck" "bag" "bcoa" "bjs" "baci" "bif" "ble" "bmet" "gst" "bacw" "bacp" "bacb" "baco" "bacy" "bacl" "balm" "beo" "bsm" "bsj" "bon" "bgy" "bfx" "bgi" "bwh" "bxi" "bhk" "bkw" "bbev" "bko" "balt" "bacs" "bmur" "bsaf" "bit" "bacq" "oih" "ocn" "gka" "gte" "gtk" "gtm" "gli" "gtn" "gwc" "gyc" "gya" "gct" "gmc" "ggh" "gjf" "gea" "gel" "gse" "gsr" "gej" "gth" "ptl" "afl" "agn" "anm" "aamy" "anl" "axl" "lsp" "lgy" "lfu" "lys" "lyb" "lyz" "lyg" "hhd" "hmn" "hli" "tap" "vir" "vhl" "vig" "vil" "vne" "vpn" "lao" "fpn" "far" 
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