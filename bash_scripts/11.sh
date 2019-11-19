 
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
  "xoy" "xoz" "xpo" "ypf" "ypg" "yph" "ypi" "ypj" "ypk" "ypl" "ypm" "ypn" "ypo" "ypp" "ypq" "ypr" "yps" "ypt" "ypu" "ypv" "ypw" "ypx" "ypy" "ypz" "yrb" "bceo" "bcep" "bcew" "bch" "bchi" "bci" "bcib" "bcig" "bcj" "bck" "bcl" "bcm" "bcn" "bcoa" "bcom" "bcon" "bprl" "bpsa" "bpsd" "bpse" "bpsh" "bpsi" "bpsl" "bpsm" "bpso" "bpsu" "bpt" "bpu" "bpum" "bpus" "bpw" "bpx" "bpy" "bpz" "bql" "bqr" "bqu" "btrh" "btrm" "bts" "btv" "btx" "bty" "btz" "bua" "bub" "buc" "bud" "bue" "buf" "bug" "bui" "buk" "buo" "bup" "buq" "bur" "cdn" "cdq" "cdu" "cea" "ced" "ceh" "cek" "cel" "cell" "cpsd" "dax" "dfa" "dnx" "ecan" "eclz" "ecz" "eta" "frc" "gps" "gtn" "hpx" "bmae" "pnr" "lez" "lfa" "lfe" "lff" "lfu" "lga" "lgi" "lgl" "lgn" "lgt" "lgu" "lgy" "lha" "lhd" "lhe" "lhf" "lhh" 
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