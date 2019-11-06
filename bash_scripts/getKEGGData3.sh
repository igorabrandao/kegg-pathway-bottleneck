 
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
  "aor" "ang" "afv" "act" "nfi" "pcs" "pdp" "cim" "cpw" "pbl" "pbn" "ure" "abe" "tve" "aje" "pno" "pte" "bze" "bsc" "bor" "aalt" "ztr" "pfj" "bcom" "npa" "tml" "spo" "cne" "cnb" "cgi" "tms" "ppl" "tvs" "dsq" "pco" "shs" "hir" "psq" "adl" "fme" "gtr" "lbc" "mpr" "mrr" "cci" "scm" "abp" "abv" "cput" "sla" "wse" "wic" "uma" "pfp" "mgl" "mrt" "pgr" "mlr" "ecu" "ein" "ehe" "ero" "nce" "mbr" "sre" "ddi" "dpp" "dfa" "ehi" "edi" "eiv" "acan" "pfa" "pfd" "pfh" "pyo" "pcb" "pbe" "pkn" "pvx" "pcy" "tan" "tpv" "tot" "beq" "bbo" "bmic" "cpv" "cho" "tgo" "tet" "ptm" "smin" "pti" "fcy" "tps" "aaf" "ngd" "pif" "psoj" "spar"
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