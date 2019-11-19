 
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
  "yin" "ykr" "yma" "ypa" "ypac" "ypb" "ypc" "ypd" "reb" "rec" "red" "rei" "ren" "rep" "ret" "rez" "rfe" "rfr" "rga" "rge" "rgl" "rhc" "rhd" "rhe" "rhf" "act" "actn" "actt" "acu" "acum" "acv" "acx" "acz" "ade" "adf" "adh" "adl" "adt" "anx" "aoa" "aof" "aoh" "aon" "aor" "apa" "apac" "apag" "apb" "apc" "apd" "ape" "apf" "apg" "aph" "apha" "api" "apj" "apk" "apl" "apm" "apor" "bmaf" "bmai" "bmal" "bman" "bmaq" "bmay" "bmaz" "bmb" "bmc" "bme" "bmec" "bmee" "bmel" "bmet" "bmf" "bmg" "bmi" "bmj" "bmk" "cir" "cit" "cja" "cjb" "cjd" "cjej" "cjen" "cjer" "cjeu" "cji" "cjj" "cjl" "cjm" "cjn" "cjp" "cjq" "cjs" "cjv" "cjw" "cjx" "cjy" "cjz" "cko" "clec" "clg" "cpf" "cpj" "cpm" "cpot" "cpr" "cpro" "cprv" "cps" "dfg" "dfl" "dfn" "dgg" "dgt" "dha" "dhk" "dhy" "dia" 
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