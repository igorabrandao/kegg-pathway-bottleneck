 
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
  "sig" "sfd" "sip" "spf" "bby" "oca" "sug" "exm" "azm" "ssan" "stv" "spa" "bwd" "suf" "exu" "azz" "snj" "spat" "spb" "bcl" "six" "rsk" "same" "rcp" "stg" "saua" "stra" "got" "smy" "tmo" "oco" "sino" "bpu" "rhp" "stx" "strn" "gmo" "saue" "span" "thal" "bop" "bpum" "rbl" "soz" "ssob" "geq" "saun" "bpus" "jan" "stz" "srq" "bbe" "saus" "efk" "rde" "seqi" "blr" "spyh" "sauu" "bos" "eah" "lpl" "bfm" "spyo" "saug" "skr" "bvv" "esj" "txi" "lpj" "bagr" "spn" "sauz" "splm" "bmq" "rli" "saut" "pjd" "thac" "atu" "boi" "splk" "lps" "vgo" "sauk" "lpz" "ppy" "spkc" "magq" "atf" "bhe" "pcon" "sauv" "ppo" "bqu" "sphf" "spx" "avi" "ljf" "bck" "pzh" "sauw" "bqr" "ppol" "snt" "ncb" "spha" "paro" "saux" "bbk" "ppq" "snd" "fer" "ljh" "agr" "btr" "ppoy" "bag" "spnn" "pbr" "ljn" "agc" "btx" "pms" "bcoa" "sne" "mgm" "sauy" "lac" "aro" "paru" "spv" "pub" "sauf" "sjp" "lai" "agt" "pamn" "sab" "pel" "sch" "lad" "ret" "pmut" "bjs" "pmq" "apc" "suy" "ssy" "laf" "rec" "pars" "baci" "pmw" "apm" "saub" "bcd" "syb" "snm" "lsa" "parr" "bif" "pta" "dsh" "ble" "plv" "baus" "kvu" "bmet" "spp" "lsl" "ecog" "psab" "saum" "rep" "lsi" "mai" "pdu" "sauc" "bvn" "spmi" "lsj" "man" "sni" "saur" "kvl" "gst" "banc" "spng" "kro" "bacw" "bapi" "snb" "bacp" "psf" "ldb" "bart" "pgv" "sphb" "saui" "rle" "pgm" "lbu" "bara" "snp" "bacb" "phr" "pga" "barw" "snx" "sphr" "saud" "pod" "rlt" "barr" "snu" "sams" "sinb" "paen" "rlg" "baro" "spne" "ldl" "suh" "spht" "paef" "baco" "pgl" "pstg" "sep" "shyd" "paeq" "bacy" "pgd" "apb" "barj" "lbr" "spnu" "ser" "sya" "pste" "bacl" "rlb" "lbk" "sclo" "paea" "balm" "sepp" "rlu" "lca" "xau" "spno" "spyg" "paee" "beo" "php" "rtr" "afe" "ppic" "rir" "afr" "sag" "phq" "rhl" "acu" "sha" "paeh" "cij" "bsm" "azc" "san" "bsj" "acz" "sno" "shh" "lpq" "oat" "sak" "bon" "afi" "lne" "lpi" "ssp" "oar" "sphg" "paej" "afj" "bgy" "sca" "lpap" "otm" "sfla" "pbj" "bfx" "maes" "mea" "slg" "lcb" "sgc" "oct" "sphy" "pih" "rhn" "lcs" "sln" "lmd" "blas" "pri" "rpha" "lce" "ssd" "lej" "ppeo" "mfn" "sags" "bgi" "mdi" "htl" "bwh" "sagm" "sdt" "lcw" "rht" "laqu" "pnp" "rdi" "lcl" "red" "pow" "rhx" "smic" "sagi" "bxi" "mpo" "swa" "bsr" "ptp" "pbv" "rhv" "sphs" "sagr" "bhk" "mza" "spas" "pxl" "rhk" "eli" "lcx" "sagp" "pyg" "rez" "elq" "lga" "sagc" "pswu" "rjg" "ery" "lre" "sagt" "sxy" "bkw" "cmar" "mrd" "bsh" "ceh" "bsy" "lrf" "sage" "pdh" "sxl" "egn" "rhr" "cmag" "mno" "lru" "sagg" "pib" "sxo" "efv" "bbev" "ngl" "sagn" "pcx" "shu" "erk" "bko" "ngg" "smu" "pkb" "scap" "err" "balt" "lrt" "neo" "malg" "bsul" "mor" "meta" "smc" "aay" "lrr" "nen" "bacs" "paih" "smut" "amx" "sscz" "lhe" "las" "bmur" "pvo" "bss" "rsu" "maqu" "sagq" "lhl" "laa" "bsaf" "plw" "bst" "rhm" 
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