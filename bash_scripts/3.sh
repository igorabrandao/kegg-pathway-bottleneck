 
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
  "pacr" "xom" "asg" "tni" "hia" "spsy" "lpf" "ten" "selo" "fta" "pbar" "pgm" "pxl" "vmo" "bmr" "tbw" "elp" "ahu" "lra" "bcas" "enc" "psec" "avl" "lba" "ypg" "suld" "sot" "cbm" "pet" "yep" "bci" "aaf" "rri" "actt" "hei" "barw" "btra" "kpa" "lmoq" "msl" "asj" "plum" "hee" "spyg" "soz" "ecq" "adh" "paca" "nit" "bao" "pman" "psh" "leo" "gbi" "poi" "spe" "cgw" "mphy" "shew" "ebs" "bacl" "ypr" "aoa" "sla" "ptx" "bok" "frm" "bstl" "pbm" "mil" "hhf" "tat" "gtk" "hht" "sind" "bsuc" "pmib" "nco" "vra" "gti" "fhw" "bpso" "nmw" "cay" "tan" "peu" "pvo" "yti" "erk" "pstw" "lmr" "ehl" "ahs" "prq" "hbi" "aid" "cri" "same" "aue" "vch" "lee" "aec" "lcl" "tar" "dar" "lue" "thes" "btx" "spng" "lbl" "woo" "yen" "sge" "deq" "smin" "pds" "zpl" "gbc" "ecz" "eol" "rmi" "nmu" "rtn" "lca" "xbv" "slt" "vbs" "mbg" "bpq" "efk" "pll" "haa" "pbl" "sni" "tmo" "bsed" "lpil" "coh" "ecol" "syl" "col" "bho" "hdt" "paro" "hms" "aacn" "pbn" "rra" "phr" "sega" "apa" "elu" "slim" "skn" "spor" "tms" "gur" "rbg" "mev" "ecoi" "ani" "bex" "fpp" "hse" "bxh" "rpe" "eoh" "lam" "cpf" "fil" "rpt" "hco" "sne" "btha" "kuy" "jaz" "spas" "bfz" "vvi" "rid" "pab" "hyo" "pte" "nsi" "ili" "cox" "tot" "lef" "stx" "pary" "sino" "sfg" "cacn" "yan" "cpaf" "scl" "dal" "nag" "ppic" "paei" "nok" "gsk" "hyd" "bpse" "rmd" "ebe" "pcd" "micr" "rrz" "bajc" "hap" "amx" "fpal" "pjd" "mea" "eso" "sgo" "rsp" "cic" "koc" "bamn" "bmw" "rpq" "vca" "vbo" "elc" "hcd" "hcn" "rhh" "bcig" "aio" "ypy" "pri" "nmx" "lio" "hcs" "bara" "lmog" "dax" "pfe" "hef" "rlt" "mge" "vcl" "bmab" "cpst" "sagi" "hih" "phz" "paln" 
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