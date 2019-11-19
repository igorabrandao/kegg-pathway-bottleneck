 
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
  "mphy" "lhr" "lat" "bit" "plen" "smj" "aep" "bso" "rhc" "mee" "anh" "pphr" "metd" "lhv" "smua" "tco" "ssif" "lhh" "bacq" "stc" "anx" "scv" "lcc" "lhd" "oih" "hat" "metx" "bsq" "stl" "alb" "lar" "lfe" "ocn" "lap" "mets" "bsx" "ste" "cna" "lau" "gka" "lagg" "asoc" "meti" "bsp" "gte" "slz" "labr" "coh" "bli" "stn" "gtk" "scoh" "labp" "cohn" "miv" "lff" "bld" "cman" "daa" "saca" "lrh" "blh" "pns" "stu" "yan" "msl" "lrg" "bay" "porl" "sscu" "gtm" "abaw" "gli" "bme" "stw" "gtn" "bmel" "lrl" "snl" "phz" "baq" "mtun" "sthe" "gwc" "bmi" "sths" "gyc" "bmz" "lra" "skl" "bts" "pot" "suam" "bya" "mlg" "ssa" "gya" "bbar" "lro" "bmg" "sfq" "ssb" "gct" "chel" "shom" "bmw" "ssu" "kyr" "bamp" "spse" "gox" "bmee" "smus" "ssv" "tum" "sulz" "goh" "mcl" "bmf" "ssi" "lcr" "cdq" "suli" "goy" "gmc" "mcak" "bmb" "sss" "lam" "hdn" "suld" "gal" "ggh" "macr" "bmc" "ssf" "tab" "lay" "hdt" "bama" "ssw" "siv" "lbh" "hmc" "bamn" "sup" "ssil" "shv" "hni" "gti" "bamb" "gjf" "baa" "sdo" "sbac" "rva" "gbe" "bamt" "gea" "babo" "don" "lmo" "phl" "gbh" "bamy" "babr" "lke" "sob" "bmp" "babt" "lrm" "pln" "bao" "babb" "lsn" "lmn" "gse" "gbc" "fil" "gsr" "gbs" "fiy" "law" "gej" "acr" "deq" "pku" "ssuy" "babu" "lmy" "tom" "bql" "lho" "ssk" "babs" "lmt" "bxh" "amv" "gth" "dei" "prt" "gdi" "ptl" "lmu" "dea" "ssq" "bami" "pll" "lmoc" "babc" "thw" "bms" "lae" "sui" "afl" "gdj" "pana" "bsi" "bamc" "bvr" "lmoe" "lgn" "suo" "agn" "gxy" "pdg" "bamf" "blag" "lmob" "rmm" "lmod" "rok" "lmow" "srp" "lko" "anm" "phc" "gxl" "rhz" "bae" "pmar" "kna" "mmyr" "ssut" "lmoq" "rid" "ppla" "aamy" "keu" "yti" "ssui" "lmr" "lvs" "lku" "pfae" "anl" "ksc" "msc" "bsup" "bvm" "axl" "apt" "mbry" "bsuv" "bha" "lsp" "apw" "mros" "lgl" "bsuc" "sgo" "ban" "lmom" "plx" "aht" 
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