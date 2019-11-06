 
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
  "ppt" "ppb" "ppi" "ppx" "ppuh" "pput" "ppun" "ppud" "pfv" "pmon" "pmot" "pmos" "ppj" "por" "pst" "psb" "psyr" "psp" "pamg" "pci" "pavl" "pfl" "pprc" "ppro" "pfo" "pfs" "pfe" "pfc" "pfn" "ppz" "pfb" "pman" "ptv" "pcg" "pvr" "pazo" "poi" "pfw" "pff" "pfx" "pen" "psa" "psz" "psr" "psc" "psj" "psh" "pstu" "pstt" "pbm" "pba" "pbc" "ppuu" "pdr" "psv" "psk" "pkc" "pch" "pcz" "pcp" "pfz" "plq" "palk" "prh" "psw" "ppv" "pses" "psem" "psec" "ppsy" "psos" "pkr" "pfk" "panr" "ppsl" "pset" "psil" "pym" "pade" "psed" "avn" "avl" "avd" "acx" "pbb" "par" "pcr" "prw" "pso" "pur" "pali" "pspg" "psyg" "psyc"
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