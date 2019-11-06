 
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
  "dsi" "dya" "dan" "dsr" "dpo" "dpe" "dmn" "dwi" "dgr" "dmo" "daz" "dnv" "dhe" "dvi" "mde" "lcq" "aga" "aag" "aalb" "cqu" "ame" "bim" "bter" "ccal" "obb" "soc" "mpha" "aec" "acep" "pbar" "vem" "hst" "dqu" "cfo" "lhu" "pgc" "obo" "pcf" "nvi" "csol" "mdl" "tca" "dpa" "atd" "nvl" "bmor" "bman" "dpl" "pmac" "prap" "haw" "tnl" "pxy" "api" "dnx" "ags" "rmd" "btab" "clec" "phu" "zne" "fcd" "dpx" "pvm" "isc" "tut" "dpte" "cscu" "ptep" "cel" "cbr" "bmy" "loa" "nai" "tsp" "hro" "lgi" "pcan" "crg" "myi" "obi" "lak" "smm" "shx" "ovi" "egl" "nve" "epa" "adf" "amil" "pdam" "spis" "dgt" "hmg" "tad" "aqu" "ath" "aly" "crb" "csat" "eus" "brp" "bna" "boe" "rsz" "thj" "cpap" "cit" "cic" "tcc" "gra" "ghi" "gab" "dzi" "egr" "gmx" "gsj" "pvu" "vra" "var" "vun" "ccaj"
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