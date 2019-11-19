 
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
  "bmc" "dosa" "llo" "aur" "geb" "ott" "ypf" "atr" "lmn" "pspi" "tmc" "lcp" "sfk" "amah" "mesc" "hcb" "xne" "synw" "ase" "masy" "kgy" "pphe" "cdm" "bmp" "gpi" "aap" "boz" "hik" "sbb" "sprn" "sbt" "sox" "psom" "pga" "gdi" "aalt" "ein" "mhaz" "kci" "eam" "bamc" "lha" "rbd" "syc" "lau" "bmai" "snt" "tal" "bbk" "sbl" "pars" "abp" "hsd" "ini" "fmi" "asu" "tau" "vhr" "tca" "vvy" "bthe" "cmed" "bgi" "mpha" "lcw" "mmyr" "ssal" "bba" "suw" "cme" "ecog" "rrm" "palk" "saum" "oat" "senl" "rpy" "xau" "mcb" "sdy" "vei" "ypj" "ssut" "scv" "sbag" "scla" "cfi" "pob" "tmar" "age" "elr" "paih" "ftt" "pfd" "rpk" "eum" "cot" "bgm" "sly" "bss" "sat" "ocn" "aaqu" "nfa" "bif" "pku" "mmk" "yma" "ngk" "rpha" "spop" "cov" "pcan" "rpla" "micc" "tsn" "nwe" "laq" "vfl" "fgl" "lne" "care" "hyf" "ypz" "kva" "thi" "amar" "sfr" "bapf" "foh" "bub" "rac" "mey" "sfw" "psem" "bae" "sul" "vcy" "cfus" "sel" "nei" "scu" "nam" "atw" "rmc" "cjy" "gte" "dod" "sck" "slg" "aalg" "udi" "coz" "ypm" "echw" "elq" "vem" "lsd" "thd" "bcon" "app" "ssif" "dzi" "clg" "asoc" "bsc" "ngg" "bsw" "hpf" "fki" "erw" "ssa" "sphm" "pfp" "ecl" "obr" "lor" "bart" "eta" "acin" "bamy" "baab" "nau" "shx" "pus" "smj" "ljn" "peq" "kur" "cfar" "pti" "tam" "blb" "blr" "ark" "lre" "arj" "camz" "act" "nvr" "anl" "vfu" "phl" "rht" "mbd" "pwo" "vow" "pspw" "plv" "ssz" "smu" "bty" "sauc" "btrm" "tpn" "lcu" "mcj" "sphs" "lua" "ccoc" "ptq" "cspf" "labp" "cif" "thal" "uma" "pef" "oca" "shw" "sagt" "bve" "bkw" "fer" "pbt" "balt" "pyo" "acom" "rop" "rse" "dja" "alc" "spa" "pur" "spnu" "ahl" "var" "bmaz" "rmb" "pseg" "bgl" "smus" "ksc" "hak" "soi" "nan" "dvu" "bmx" "bol" "hpyk" "dav" "lyg" "rpx" "labr" "wwe" "gct" "emr" "jre" "pnc" "lal" "seo" "hes" "pdi" "rpo" "hahe" "lgu" 
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