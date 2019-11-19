 
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
  "aao" "bsj" "psi" "sek" "spnn" "gme" "eclc" "plu" "abe" "rbz" "ble" "vco" "bpw" "mmaa" "ebc" "ppuh" "rok" "tas" "cso" "shl" "sra" "hat" "eld" "pms" "nmm" "sbw" "aah" "blo" "baft" "lpi" "jan" "scou" "sss" "keu" "saue" "phei" "psp" "tee" "bsul" "nin" "pyc" "paef" "lpe" "sgc" "pset" "rbe" "mlo" "mali" "yal" "vpa" "enl" "bant" "bafe" "bak" "kzo" "sku" "bai" "dct" "buc" "thap" "rcc" "tea" "thin" "bvv" "vag" "echj" "pcon" "xfl" "mhal" "psq" "rlb" "pep" "raf" "rpw" "wri" "hpas" "hhu" "kpz" "hcc" "lsv" "sll" "epa" "lab" "cjd" "xor" "cjb" "bcoa" "bsi" "swp" "hel" "gni" "pphr" "mro" "achr" "pap" "leh" "aqu" "hst" "yfr" "abo" "shf" "cgr" "sbac" "nbg" "eck" "srb" "scap" "lul" "rup" "pand" "lsj" "hpys" "chel" "axy" "lej" "dhy" "psc" "bmul" "spf" "bso" "rpd" "mct" "pbc" "ano" "seh" "laci" "eae" "pso" "kpo" "amp" "pkc" "pah" "vcj" "agl" "amv" "ege" "bami" "vaf" "bpc" "sdl" "bmee" "neu" "lmob" "aof" "ssd" "bpsi" "samy" "ceh" "rpv" "hpyo" "smul" "lbh" "rhm" "pvx" "aes" "hhh" "pstu" "nen" "hct" "ecw" "tgr" "tig" "aem" "lip" "psf" "mbw" "kpr" "mmae" "acp" "bpsd" "bpa" "yin" "sphr" "mee" "wgl" "rhe" "span" "sauy" "bsk" "mla" "bprl" "tco" "hhz" "lsa" "kie" "cqi" "ype" "sob" "ppoa" "ppoy" "bmec" "sbi" "ecy" "fam" "exm" "sec" "bapu" "part" "cjs" "rtb" "pvm" "lmod" "hpr" "smut" "lff" "halz" "crt" "mip" "lrg" "sbf" "bih" "pdu" "senn" "slz" "lgt" "aln" "kvl" "ram" "lcr" "egl" "nwr" "hbe" "cpsd" "clu" "sezo" "tmm" "opf" "mhi" "magq" "apha" "pana" "csol" "bbs" "alt" "ecos" "mpur" "pfh" "cfo" "cbg" "ria" "ngo" "vtu" "ecla" "mrd" "pagc" "apro" "ypq" "mcat" "cthr" "acra" "psu" "hix" "dte" "npl" "taj" "van" "odi" "lps" "kaf" "pcg" "ear" "psd" "drt" "rrc" "pgo" "ero" "tab" "afa" "ntt" "xoo" "eal" "kos" "cfx" "gth" "llg" "gsl" "eah" "bdi" "sths" "lsi" "sth" "mnr" "xfs" "vna" "pul" "amim" "ftg" "kgo" "bpm" "lke" "vvu" "asw" "masz" "ske" "lho" "wso" "hvn" "ebu"
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