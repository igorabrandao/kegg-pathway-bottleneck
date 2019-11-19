 
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
  "hpt" "hpu" "sfq" "sfr" "sfu" "sfw" "sfz" "sgc" "sge" "sgl" "sgn" "sgo" "sgp" "paf" "pag" "pagc" "pagg" "pah" "paih" "paj" "palh" "palk" "paln" "paly" "pam" "pamg" "pami" "pamn" "pan" "pana" "panc" "pand" "panp" "panr" "pans" "pant" "pap" "agc" "age" "agi" "agl" "agn" "ago" "agr" "ags" "agt" "ahh" "ahi" "ahj" "ahl" "ahn" "ahs" "aht" "ahu" "ahz" "aid" "aio" "aje" "ajn" "dba" "dco" "dcr" "dct" "dda" "ddc" "ddi" "ddn" "ddq" "dds" "dea" "def" "dei" "del" "den" "deo" "deq" "des" "deu" "ppsy" "ppt" "ppuh" "ppun" "pput" "ppuu" "ppv" "ppw" "ppx" "ppy" "ppz" "prap" "prf" "mpur" "mpy" "mrd" "mrh" "mrm" "mro" "mros" "mrr" "mrs" "mrt" "msc" "msd" "mse" "msl" "hyd" "hye" "hyf" "hyo" "hys" "icp" "idt" "ili" "bon" "boo" "bop" "bor" "bos" "bov" "boz" "bpa" "bpc" 
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