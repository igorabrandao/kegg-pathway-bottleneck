 
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
  "dco" "amac" "gem" "ech" "bfx" "brb" "dfg" "smo" "ecan" "dau" "blat" "bya" "pamg" "caa" "sfu" "fad" "dun" "lmow" "bpk" "acz" "paea" "cbr" "lfe" "lsp" "spoi" "mass" "bgd" "mgm" "ape" "klw" "eiv" "dto" "gea" "alb" "echg" "hpyh" "lmoj" "baq" "cbn" "bpz" "lil" "sphb" "bacp" "aeo" "ebl" "cjw" "rhu" "saca" "tcx" "snz" "esj" "cpm" "ack" "pff" "mab" "halo" "xci" "yps" "sxy" "mdn" "mpsy" "hjt" "ecp" "bmq" "rrn" "pdag" "sbj" "mdi" "rhn" "rsc" "kpc" "nmn" "rmm" "btq" "cjm" "kpb" "nme" "hha" "ftl" "salv" "apc" "sre" "kpm" "fcy" "hpyg" "kqv" "hci" "skl" "tbn" "yph" "kll" "hol" "mer" "btei" "bqu" "tom" "cbu" "isc" "shj" "bcs" "porl" "pcp" "palh" "pmq" "mib" "sphg" "vmi" "bmor" "smer" "afl" "nob" "rjg" "cha" "llu" "dba" "eln" "sur" "pfs" "spm" "svo" "mej" "wen" "hce" "apm" "ypd" "rdi" "rfr" "alp" "afi" "balm" "mema" "arb" "kpnk" "psab" "prj" "mthe" "tmi" "asy" "bcl" "xcc" "pcy" "ngd" "saln" "haz" "fek" "saga" "cvc" "meh" "mse" "kok" "lbr" "lph" "amed" "del" "xho" "apt" "let" "mpru" "pbv" "lcx" "xdo" "ced" "xoz" "hax" "dak" "bha" "amk" "bbh" "pca" "bthy" "tcc" "senr" "bmg" "pbor" "pdh" "xfu" "ppv" "see" "kph" "eih" "spno" "wpa" "spos" "splk" "lcc" "bbo" "goh" "ypu" "axe" "bbat" "blp" "edj" "bas" "nmp" "nce" "mem" "dpi" "ftr" "gal" "yef" "kox" "sele" "ssob" "bph" "spht" "lbu" "afq" "swa" "saua" "laa" "efe" "rde" "bmu" "ocd" "pgn" "sev" "hro" "ldl" "selt" "bfm" "acu" "fsi" "ecor" "ppa" "mbry" "cui" "siv" "adv" "geh" "ppun" "dao" "babr" "sua" "stw" "mvs" "cavi" "ypx" "ppsy" "srs" "dsq" "dmi" "snx" "nmh" "gao" "buk" "fnf" "ehe" "ddq" "mhat" "btv" "kma" 
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