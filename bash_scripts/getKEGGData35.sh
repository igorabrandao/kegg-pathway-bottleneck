 
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
  "mtc" "mra" "mtf" "mtb" "mtk" "mtz" "mtg" "mti" "mte" "mtur" "mtl" "mto" "mtd" "mtn" "mtj" "mtub" "mtuc" "mtue" "mtx" "mtuh" "mtul" "mtut" "mtuu" "mtq" "mbo" "mbb" "mbt" "mbm" "mbk" "mbx" "maf" "mmic" "mce" "mcq" "mcv" "mcx" "mcz" "mle" "mlb" "mpa" "mao" "mavi" "mavu" "mav" "mit" "mia" "mid" "myo" "mir" "mchi" "mmal" "mlp" "msa" "mul" "mmc" "mkm" "mjl" "mmi" "mmae" "mmm" "mli" "mkn" "myv" "mye" "mhad" "myn" "mdx" "mshg" "msm" "msg" "msb" "msn" "msh" "mva" "mgi" "msp" "mcb" "mne" "mgo" "mft" "mphl" "mvq" "mll" "mrh" "mthn" "mhas" "mab" "mmv" "mabb" "mabl" "mche" "miz" "mste" "msao" "msal" "mjd" "mter" "asd" "cgl" "cgb" "cgu" "cgt" "cgs" "cgg" "cgm" "cgj" "cgq" "cgx" "cef" "cdi" "cdp" "cdh" "cdt" "cde" "cdr" "cda" "cdz" "cdb" "cds" "cdd" "cdw" "cdv" "cdip" "cjk" 
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