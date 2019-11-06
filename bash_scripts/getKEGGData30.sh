 
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
  "sda" "sdc" "sdq" "sga" "sgg" "sgt" "smb" "sor" "stk" "stb" "scp" "scf" "ssr" "stf" "stj" "strs" "ssah" "std" "smn" "sif" "sie" "sib" "siu" "sang" "sanc" "sans" "scg" "scon" "scos" "soi" "sik" "siq" "sio" "siz" "slu" "sig" "sip" "stv" "spat" "stra" "strn" "ssob" "srq" "seqi" "lpl" "lpj" "lpt" "lps" "lpr" "lpz" "lpb" "ljo" "ljf" "ljh" "ljn" "lac" "lai" "lad" "laf" "lsa" "lsl" "lsi" "lsj" "ldb" "lbu" "lde" "ldl" "lbr" "lbk" "lca" "lcz" "lpq" "lpi" "lpap" "lcb" "lcs" "lce" "lcw" "lcl" "lcx" "lga" "lre" "lrf" "lru" "lrt" "lrr" "lhe" "lhl" "lhr" "lhv" "lhh" "lhd" "lfe" "lfr" "lff" "lrh" "lrg" "lrl" "lra" "lro" "lrc" "lcr" "lam" "lay" "lbh" "lbn" "lke" "lrm" "lsn" "law" "lho" "lmu" "lae" "lgn" "lko" "lhi" "lku" "lgl" "lpx" "lor" "lpar" "lle" "lpd" "lje" "lcu" "lct" "lbt" 
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