 
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
  "cur" "cua" "car" "ckp" "cpu" "cpl" "cpg" "cpp" "cpk" "cpq" "cpx" "cpz" "cor" "cop" "cod" "cos" "coi" "coe" "cou" "cpse" "cpsu" "cpsf" "crd" "cul" "cuc" "cue" "cun" "cus" "cuq" "cuz" "cuj" "cva" "chn" "ccn" "cter" "cmd" "caz" "cfn" "ccg" "cvt" "cgy" "cax" "cii" "cuv" "coa" "cdo" "chm" "csx" "cmq" "cku" "ccj" "cmv" "cei" "cted" "cut" "clw" "cdx" "csp" "csta" "ccjz" "cfk" "cpho" "cfc" "cgv" "cstr" "caqu" "csph" "camg" "cmin" "cpeg" "bfv" "nfa" "nfr" "ncy" "nbr" "nno" "nsl" "nsr" "ntp" "noz" "nod" "rha" "rer" "rey" "reb" "rop" "roa" "req" "rpy" "rhb" "rav" "rfa" "rhw" "rhs" "rrz" "rhu" "rqi" "rhq" "rhod" "rrt" "rby" "gbr" "gpo" "gor" "goq" "gta" "goc" "git" "gru" "gom" "gav" "tpr" "tsm" "srt" "dtm" "dit" "diz" "dpc" "dlu" "cbq" "toy" "sco" "salb" "sma" "sgr" "sgb" "scb" 
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