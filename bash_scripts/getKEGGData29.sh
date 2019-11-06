 
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
  "pgq" "ntr" "lfb" "tvu" "lla" "llk" "llt" "lls" "lld" "llx" "llc" "llm" "llr" "lln" "lli" "llw" "llj" "lgr" "lgv" "lpk" "lrn" "lact" "spy" "spz" "spym" "spya" "spm" "spg" "sps" "sph" "spi" "spj" "spk" "spf" "spa" "spb" "stg" "stx" "soz" "stz" "spyh" "spyo" "spn" "spd" "spr" "spw" "sjj" "snv" "spx" "snt" "snd" "spnn" "sne" "spv" "snc" "snm" "spp" "sni" "spng" "snb" "snp" "snx" "snu" "spne" "spnu" "spnm" "spno" "sag" "san" "sak" "sgc" "sags" "sagl" "sagm" "sagi" "sagr" "sagp" "sagc" "sagt" "sage" "sagg" "sagn" "smu" "smc" "smut" "smj" "smua" "stc" "stl" "ste" "stn" "stu" "stw" "sthe" "sths" "ssa" "ssb" "ssu" "ssv" "ssi" "sss" "ssf" "ssw" "sup" "ssus" "sst" "ssuy" "ssk" "ssq" "sui" "suo" "srp" "ssut" "ssui" "sgo" "sez" "seq" "sezo" "sequ" "seu" "sub" "sds" "sdg" 
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