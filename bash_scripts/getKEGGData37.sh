 
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
  "ssx" "svl" "sct" "scy" "sfa" "sbh" "shy" "sho" "sve" "sdv" "sals" "strp" "sfi" "sci" "src" "salu" "sall" "slv" "sgu" "svt" "stre" "scw" "sld" "slc" "sxi" "strm" "strc" "samb" "spri" "scz" "scx" "srw" "strf" "sle" "srn" "spav" "strt" "sclf" "sgs" "stsi" "sls" "snr" "splu" "strd" "snw" "sauo" "ssia" "svu" "spun" "sgv" "smal" "slau" "salf" "salj" "slx" "stro" "sfk" "snz" "sge" "ksk" "kab" "kau" "kit" "stri" "twh" "tws" "lxx" "lxy" "cmi" "cms" "cmc" "cmh" "ccap" "mts" "mim" "mio" "mix" "mip" "mcw" "mpal" "mih" "micr" "maur" "mhos" "mfol" "moo" "rla" "rpla" "rtx" "rtc" "rtn" "rry" "ria" "rfs" "cum" "cub" "cug" "mvd" "frp" "agy" "agm" "agf" "cart" "cry" "cphy" "amin" "aum" "auw" "psai" "malk" "myl" "salc" "sala" "sald" "hum" "huw" "gry" "lyd" "lyk" "plap" "leu" 
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