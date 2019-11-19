 
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
  "dic" "dja" "dji" "dsl" "dso" "dsq" "dsu" "dte" "dti" "dto" "dtr" "dts" "dtx" "dun" "dvg" "dvl" "dvu" "dye" "dzc" "dze" "dzi" "eaa" "eab" "eae" "eah" "ecg" "ech" "echg" "echj" "echl" "echs" "echv" "echw" "eci" "eck" "ecl" "ecla" "eclc" "ecle" "eclg" "eclo" "ecls" "ecly" "edh" "edi" "edj" "eec" "eex" "efe" "efk" "efv" "hak" "halo" "halz" "ham" "han" "hap" "hat" "jre" "kab" "kaf" "kak" "kbl" "kbt" "kci" "kct" "kde" "keu" "kga" "kgo" "kgy" "kie" "kit" "kla" "kle" "kll" "kma" "kmi" "kmx" "kna" "koc" "koe" "kok" "kos" "kox" "kpa" "kpb" "kpc" "kpe" "kpg" "kph" "kpi" "kpk" "kpm" "kpne" "kpnk" "kpnu" "kpo" "kpq" "kpr" "kpt" "kpw" "kpy" "kpz" "kqu" "kqv" "malg" "mali" "man" "maqu" "marc" "marf" "myi" "myv" "mza" "mzh" "nac" "naf" "nag" "nai" "naj" "nam" "nan" 
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