 
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
  "lcp" "lgu" "lez" "lem" "lmb" "lyt" "lyj" "lum" "lue" "lus" "lug" "thes" "xbc" "fau" "rhd" "rgl" "dji" "dja" "dtx" "dye" "dko" "lrz" "lpy" "xba" "rbd" "vch" "vcf" "vcs" "vce" "vcq" "vcj" "vci" "vco" "vcr" "vcm" "vcl" "vcx" "vcz" "vvu" "vvy" "vvm" "vvl" "vpa" "vpb" "vpk" "vpf" "vph" "vha" "vca" "vag" "vex" "vdb" "vhr" "vna" "vow" "vro" "vsp" "vej" "vfu" "vni" "van" "lag" "vau" "vcy" "vct" "vtu" "vfl" "vmi" "vbr" "vsc" "vga" "vsh" "vqi" "vta" "vaf" "vnl" "vfi" "vfm" "vsa" "awd" "ppr" "pgb" "pds" "gho" "pmai" "saly" "sks" "pae" "paev" "paei" "pau" "pap" "pag" "paf" "pnc"
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