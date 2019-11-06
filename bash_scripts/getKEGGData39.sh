 
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
  "mgg" "tfu" "nda" "nal" "ngv" "strr" "tcu" "actw" "sro" "noa" "fra" "fre" "fri" "fal" "fsy" "ace" "nml" "nak" "gob" "bsd" "mmar" "kra" "sen" "svi" "sacc" "amd" "amn" "amm" "amz" "aoi" "aja" "amq" "amyc" "amyb" "aab" "pdx" "psea" "psee" "pseh" "pseq" "pecq" "phh" "paut" "ami" "apre" "sesp" "kal" "kphy" "led" "ahm" "acti" "acad" "ahg" "acta" "alo" "pmad" "stp" "saq" "mau" "mil" "micb" "mtua" "mich" "vma" "ams" "ase" "actn" "afs" "acts" "plk" "plab" "plat" "cai" "sna" "sale" "ahe" "mcu" "tpy" "tpyo" "tbw" "asg" "actt" "amy" "soo" "acq" "aos" "ard" "avu" "actp" "actc" "acto" "ane" "ahw" "fsl" "flh" "blo" "blj" "bln" "blon" "blf" "bll" "blb" "blm" "blk" "blg" "blz" "blx" "bad" "badl" "bado" "bla" "blc" "blt" "bbb" "bbc" "bnm" "blv" "blw" "bls" "bani" "banl" "bni" "banm" 
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