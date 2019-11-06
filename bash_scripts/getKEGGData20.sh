 
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
  "sur" "age" "mbd" "cfus" "vin" "scl" "scu" "ccro" "samy" "llu" "mrm" "hoh" "sat" "dao" "dti" "sfu" "dax" "dbr" "hmr" "dav" "bsed" "bba" "bbat" "bbw" "bbac" "bex" "bdq" "bmx" "hax" "bsto" "sbf" "rpr" "rpo" "rpw" "rpz" "rpg" "rps" "rpv" "rpq" "rpl" "rpn" "rty" "rtt" "rtb" "rcm" "rcc" "rbe" "rbo" "rco" "rfe" "rak" "rri" "rrj" "rra" "rrc" "rrh" "rrb" "rrn" "rrp" "rrm" "rrr" "rms" "rmi" "rpk" "raf" "rhe" "rja" "rsv" "rsw" "rph" "rau" "rmo" "rpp" "rre" "ram" "rab" "rmc" "ric" "ots" "ott" "ptc" "wol" "wri" "wen" "wed" "wpi" "wbm" "woo" "wcl" "weo" "wpp" "ama" "amf" "amw" "amp" "acn" "aph" "apy" "apd" "apha" "aoh" "eru" "erw" "erg" "ecn" "ech" "echa" "echj" "echl" "echs" "echv" "echw" "echp" "emr" "ehh" "nse" "nri" "nhm" "mmn" "fso" "rbt" "ren" "paca" "caq" "naf" "eaa" "mlo" 
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