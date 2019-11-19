 
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
  "hrb" "hro" "hru" "ngd" "nge" "ngg" "ngk" "ngl" "ngo" "nha" "nhl" "nhm" "nia" "nid" "nii" "nin" "nis" "nit" "nkr" "nlc" "nma" "nmc" "nme" "nmg" "nmh" "pns" "pob" "pod" "poi" "pol" "pom" "pop" "por" "porl" "pos" "pot" "pow" "ppa" "ppb" "ppd" "rsz" "rtb" "rtc" "rtn" "rtr" "rtt" "rty" "rua" "rue" "run" "rup" "rva" "rvi" "saal" "sab" "saca" "saci" "psh" "psi" "psil" "psj" "psk" "pso" "psom" "psos" "psp" "pspg" "pspi" "pspo" "pspw" "psq" "psr" "pst" "psta" "pste" "pstg" "psts" "pstu" "pstw" "psu" "hcs" "hct" "hdn" "hdt" "hdu" "heb" "hee" "hef" "hei" "hel" "hem" "hen" "hep" "erc" "erg" "erj" "erk" "ern" "ero" "err" "eru" "erw" "ery" "esa" "esc" "ese" "esi" "esj" "esl" "eso" "ilo" "jeu" "klw" "kpj" "labp" "lch" "lew" "lpo" "mbd" "meu" "mhaq" "doa" "dod" "don" 
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