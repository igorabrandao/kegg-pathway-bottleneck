 
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
  "cbx" "oto" "otk" "lim" "lih" "hyr" "hyb" "hyl" "hyc" "hpse" "dpy" "simp" "melm" "mela" "cbaa" "cbab" "mpt" "metp" "har" "mms" "jag" "jab" "jaz" "jal" "jsv" "jaj" "hse" "hsz" "hht" "hrb" "hee" "hhf" "zin" "cfu" "care" "cpra" "mnr" "masw" "mass" "masz" "mtim" "masy" "mali" "mum" "ofo" "upv" "nok" "sutt" "sutk" "lch" "tin" "thi" "rge" "rbn" "rdp" "pkt" "miu" "rgu" "aon" "snn" "bbag" "bbay" "pbh" "neu" "net" "nit" "nii" "nco" "nur" "nst" "nmu" "nlc" "shd" "metr" "tbd" "mfa" "mmb" "meh" "mei" "mep" "mbac" "mbat" "meu" "slt" "gca" "sdr" "sulf" "fam" "eba" "dsu" "rbu" "otr" "rbh" "dar" "dey" "azo" "aza" "azi" "aoa" "atw" "acom" "azd" "azr" "tmz" "thu" "tcl" "thk" "tak" "zpa" "app" "tpn" "tpq" "tpj" "kci" "kct" "kbl" "kbt" "kde" "kga" "kon" "kso" "ssdc"
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