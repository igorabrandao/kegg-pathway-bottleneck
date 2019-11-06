 
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
  "oan" "oah" "ops" "och" "bja" "bju" "bjp" "bra" "bbt" "brs" "aol" "brc" "brad" "bic" "bro" "brk" "bot" "brq" "bgq" "bgz" "rpa" "rpb" "rpc" "rpd" "rpe" "rpt" "rpx" "nwi" "nha" "oca" "ocg" "oco" "bop" "bos" "bvv" "boi" "bof" "vgo" "bhe" "bhn" "bhs" "bqu" "bqr" "bbk" "btr" "btx" "bgr" "bcd" "baus" "bvn" "banc" "bapi" "bart" "bara" "barw" "barr" "baro" "barj" "bez" "xau" "azc" "sno" "lne" "mex" "mea" "mdi" "mch" "mpo" "mza" "mrd" "met" "mno" "mor" "meta" "maqu" "mphy" "mee" "metd" "metx" "mets" "meti" "moc" "miv" "bid" "msl" "mtun" "mlg" "bbar" "chel" "cdq" "hdn" "hdt" "hmc" "hni" "rva" "phl" "fil" "fiy" "deq" "dei" "dea" "bvr" "blag" "rhz" "mmyr" "yti" "msc" "mbry" "mros" "mtw" "pleo" "mey" "maad" "mmed" "aua" "brn" "hci" "hct" "hcc" "hcd" "mcg" "thd" "psin" 
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