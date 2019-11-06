 
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
  "mtr" "cam" "lja" "adu" "aip" "lang" "fve" "rcn" "pper" "pmum" "pavi" "mdm" "pxb" "zju" "csv" "cmo" "mcha" "cmax" "cmos" "cpep" "rcu" "jcu" "hbr" "mesc" "pop" "peu" "jre" "qsu" "vvi" "sly" "spen" "sot" "cann" "nta" "nsy" "nto" "nau" "ini" "sind" "oeu" "han" "lsv" "ccav" "dcr" "bvg" "soe" "cqi" "nnu" "psom" "osa" "dosa" "obr" "bdi" "ats" "sbi" "zma" "sita" "pda" "egu" "mus" "dct" "peq" "aof" "atr" "smo" "ppp" "cre" "vcn" "mng" "olu" "ota" "bpg" "mis" "mpp" "csl" "cvr" "apro" "cme" "gsl" "ccp" "sce" "ago" "erc" "kla" "kmx" "lth" "vpo" "zro" "cgr" "ncs" "ndi" "tpf" "tbl" "tdl" "kaf" "ppa" "dha" "pic" "pgu" "spaa" "lel" "cal" "ctp" "cot" "cdu" "cten" "yli" "clu" "clus" "caur" "slb" "pkz" "ncr" "nte" "smp" "pan" "ttt" "mtm" "cthr" "mgr" "tmn" "ssck" 
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