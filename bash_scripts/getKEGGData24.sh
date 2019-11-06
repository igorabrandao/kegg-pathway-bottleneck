 
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
  "zmm" "zmb" "zmi" "zmc" "zmr" "zmp" "nar" "npp" "npn" "nre" "nov" "not" "ntd" "nor" "ngf" "sal" "sphk" "sphp" "smag" "smaz" "ster" "sgi" "sphl" "sphq" "spho" "sphu" "swi" "sphd" "sphm" "stax" "sphi" "ssan" "snj" "smy" "span" "skr" "splm" "splk" "sphj" "spkc" "sphc" "sphf" "spha" "spau" "sjp" "sch" "ssy" "syb" "sbd" "spmi" "sphb" "sphr" "sinb" "spht" "shyd" "sya" "sclo" "spyg" "cij" "sphg" "sfla" "sphy" "blas" "bfw" "rdi" "smic" "sphs" "eli" "elq" "ery" "egn" "efv" "erk" "err" "aay" "amx" "aep" "anh" "ado" "alb" "cna" "cman" "pns" "porl" "phz" "pot" "gox" "goh" "goy" "gal" "gti" "gbe" "gbh" "gbc" "gbs" "acr" "amv" "gdi" "gdj" "gxy" "gxl" "kna" "keu" "ksc" "apt" "apw" "apf" "apu" "apg" "apq" "apx" "apz" "apk" "asz" "asv" "aace" "aper" "apom" "ato" "aasc" 
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