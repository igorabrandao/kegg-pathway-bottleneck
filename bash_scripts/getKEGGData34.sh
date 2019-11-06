 
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
  "mcac" "mcap" "mcar" "mcai" "mlc" "mlh" "mmo" "mhy" "mhj" "mhp" "mhn" "mhyl" "mhyo" "msy" "mso" "maa" "mal" "mat" "mco" "mho" "mhom" "mcd" "mhr" "mhh" "mhm" "mhs" "mhv" "mfr" "mfm" "mfp" "mbv" "mbh" "mbi" "mbq" "mha" "mhf" "mss" "msk" "mpf" "mput" "mhe" "mwe" "mhl" "mcy" "mhb" "mpv" "mov" "mbc" "mcr" "mcm" "mgj" "mfq" "mcan" "myt" "mds" "mgb" "mcas" "mck" "marg" "myg" "mpul" "mbov" "mpho" "mhyv" "mboh" "mclo" "uur" "upa" "upr" "uue" "hcr" "poy" "ayw" "mbp" "pml" "pal" "nzs" "psol" "pzi" "acl" "abra" "apal" "aoc" "aaxa" "mfl" "mfw" "mchc" "mlac" "ment" "msyr" "mtab" "mcol" "elj" "esx" "efr" "eml" "scr" "ssyr" "sdi" "stai" "sapi" "smir" "smia" "scq" "ssab" "satr" "seri" "stur" "sll" "skn" "scj" "shj" "sck" "sfz" "scou" "scla" "sprn" "spit" "mbj" "tbm" "tbz" "mtu" "mtv" 
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