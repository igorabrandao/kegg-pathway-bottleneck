 
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
  "lcv" "lah" "lamy" "lagl" "lzy" "lpg" "lcy" "laca" "lalw" "lali" "lpw" "lfm" "lng" "ppe" "ppen" "pce" "pdm" "paci" "pio" "efa" "efl" "efi" "efd" "efs" "efn" "efq" "ene" "efc" "efau" "efu" "efm" "eft" "ehr" "ecas" "emu" "edu" "ega" "ess" "eth" "egv" "eav" "mps" "mpx" "thl" "too" "tkr" "vte" "vpi" "vac" "vao" "ooe" "oen" "osi" "lme" "lmm" "lmk" "lci" "lki" "lec" "lcn" "lgs" "lge" "llf" "lgc" "lsu" "wko" "wce" "wct" "wci" "wcb" "wjo" "wpa" "wcf" "wso" "aur" "aun" "aui" "asan" "acg" "avs" "auh" "abae" "crn" "cml" "caw" "carc" "cdj" "marr" "jep" "jda" "jeh" "dpm" "cac" "cae" "cay" "cpe" "cpf" "cpr" "ctc" "ctet" "cno" "cbo" "cba" "cbh" "cby" "cbl" "cbk" "cbb" "cbi" "cbn" "cbt" "cbf" "cbm" "cbj" "cbe" "cbz" "cbei" "ckl" "ckr" "clj" "ccb" "cls" "clb" "csr" "cpas" "cpat" 
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