 
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
  "mln" "mci" "mop" "mam" "mamo" "meso" "mesw" "mesm" "mesp" "mes" "hoe" "aak" "amih" "pht" "rpod" "pla" "rbs" "sme" "smk" "smq" "smx" "smi" "smeg" "smel" "smer" "smd" "rhi" "sfh" "sfd" "six" "same" "sino" "ead" "eah" "esj" "atu" "ara" "atf" "ata" "avi" "agr" "agc" "aro" "agt" "ret" "rec" "rel" "rep" "rei" "rle" "rlt" "rlg" "rlb" "rlu" "rtr" "rir" "rhl" "rga" "rhn" "rpha" "rht" "rhx" "rhv" "rhk" "rez" "rjg" "rhr" "ngl" "ngg" "neo" "nen" "las" "laa" "lat" "lso" "lcc" "lar" "lau" "shz" "abaw" "bme" "bmel" "bmi" "bmz" "bmg" "bmw" "bmee" "bmf" "bmb" "bmc" "baa" "babo" "babr" "babt" "babb" "babu" "babs" "babc" "bms" "bsi" "bsf" "bsui" "bsup" "bsuv" "bsuc" "bmt" "bsz" "bsv" "bsw" "bsg" "bov" "bcs" "bsk" "bol" "bcar" "bcas" "bmr" "bpp" "bpv" "bcet" "bcee" "bvl" "bru" 
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