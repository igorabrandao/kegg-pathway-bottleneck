 
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
  "ndl" "vfg" "bprc" "beb" "beba" "hpy" "heo" "hpj" "hpa" "hps" "hhp" "hhq" "hhr" "hpg" "hpp" "hpb" "hpl" "hpc" "hca" "hpm" "hpe" "hpo" "hpi" "hpq" "hpw" "hpu" "hef" "hpf" "heq" "hex" "hpt" "hpz" "hpv" "hpx" "hen" "hph" "heg" "hpn" "hep" "heu" "hes" "hpys" "hcn" "hpd" "hey" "her" "hei" "hpya" "hpyk" "hpyo" "hpyl" "hpyb" "hpyc" "hpyd" "hpye" "hpyf" "hpyg" "hpyh" "hpyj" "hpyr" "hpyi" "hpyu" "hpym" "hem" "heb" "hez" "hhe" "hac" "hms" "hfe" "hbi" "hce" "hcm" "hcp" "hcb" "hhm" "hty" "hbl" "had" "het" "wsu" "tdn" "sua" "sku" "sulr" "cje" "cjb" "cjj" "cju" "cjn" "cji" "cjm" "cjs" "cjp" "cjej" "cjeu" "cjen" "cjei" "cjer" "cjv" "cjy" "cjq" "cjl" "cjw" "cjr" "cjd" "cjz" "cjx" "cff" "cft" "cfv" "cfx" "cfz" "camp" "cfp" "ccv" "cha" "cco" "ccoc" "cla"
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