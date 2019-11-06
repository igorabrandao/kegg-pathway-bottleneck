 
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
  "ftz" "ftm" "ftn" "ftx" "ftd" "fty" "fcf" "fcn" "fhi" "fph" "fpt" "fpi" "fpm" "fpx" "fpz" "fpj" "frt" "fna" "fnl" "frf" "fper" "fha" "frx" "frm" "frc" "fad" "fmi" "fgu" "tcx" "htr" "tcy" "tao" "thio" "mej" "mec" "cyq" "cza" "cyy" "psal" "thig" "tig" "blep" "noc" "nhl" "nwa" "nwr" "alv" "tvi" "tmb" "mpur" "tee" "ntt" "tsy" "rhh" "aeh" "hha" "hhk" "hhc" "ebs" "tgr" "tkm" "tni" "tti" "tvr" "ssal" "spiu" "sros" "aprs" "hna" "haz" "wma" "woc" "gai" "ttc" "hch" "hahe" "csa" "hel" "hcs" "hak" "ham" "hhu" "hco" "hsi" "halo" "hhh" "hbe" "hag" "haf" "halk" "hvn" "hol" "ple" "ply" "plr" "plo" "pld" "plb" "plc" "pli" "paly" "crp" "cru" "crc" "crt" "crh" "crv" "cri" "eme" "zpl" "haa" "cmai" "kus" "kma" "kuy" "paur" "abo" "adi" "apac" "aln" "axe" "kak" 
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