 
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
  "lrt" "snj" "gyc" "lct" "blj" "pos" "gtn" "ean" "pacn" "plo" "vsa" "mii" "eru" "mtun" "suam" "hgi" "lit" "pfr" "sedi" "sauk" "maqu" "thk" "dhk" "lrh" "cmp" "bvr" "lmy" "ypl" "soe" "hru" "mrs" "pagg" "mgn" "smua" "sena" "pmar" "dzc" "thj" "rhl" "hme" "arac" "nvl" "pfa" "psx" "sbs" "pvu" "dpa" "fla" "lru" "bxi" "sep" "glo" "dko" "ock" "ecr" "bos" "smp" "mdh" "tvr" "mbac" "baph" "mfn" "rsd" "asc" "csl" "dda" "shom" "icp" "yak" "oll" "vdi" "fva" "kpq" "bab" "goy" "mko" "fcd" "rqi" "set" "sut" "cfd" "bans" "setc" "asi" "sbp" "nma" "hoh" "plc" "dpte" "cpe" "heb" "hep" "hpz" "dze" "gsj" "wma" "gka" "ota" "arc" "ple" "pspo" "saf" "att" "mthr" "cyq" "kle" "pin" "oeu" "dha" "thac" "mlg" "ssq" "rms" "pfb" "cbc" "avr" "alz" "lew" "bch" "mod" "lad" "bacy" "sui" "kct" "lag" "kbt" "crv" "cgrn" "dsh" "sequ" "sita" "ath" "lby" "fso" "afd" "mvi" "mrt" "bgu" "xac" "fgo" "acr" "scj" "dgt" "aev" "ahh" "swi" "bsem" "ama" "anm" "lap" "rhx" "suf" "lum" "hdn" "vci" "rae" "bvi" "bths" "bpeu" "sln" "ern" "hfe" "laf" "sse" "bcew" "sja" "hpn" "cjer" "rsw" "bceo" "pdg" "xii" "pel" "kpi" "bana" "moi" "avi" "vcn" "pfn" "abl" "hja" "babu" "btz" "siz" "mvr" "lai" "pgv" "sdz" "bct" "hir" "prf" "nfi" "kqu" "lcs" "hpyc" "nmq" "cit" "buu" "rsu" "stl" "rcp" "saci" "lbc" "lis" "ehm" "pns" "lfa" "bld" "spyo" "hhc" "txi" "aje" "thio" "ppx" "sdr" "vcq" "ftf" "plb" "ply" "sfh" "pmx" "pzh" "barr" "ptep" "acii" "hpya" "rec" "psuw" "pub" "lsh" "hpyj" "pbw" "abal" "sros" "bpet" "acal" "pif" "stax" "baci" "ttk" "lpn" "baco" "lsl" "kab" "pgi" "hhe" "mgh" "kvu" "stem" "saur" "ypn" "fne" "sno" "cput" "saz" "nto" "sio" "ecls" "oct" "lrf" "apk" "csj" "pfae" "senq" "bamt" "plw" "bue" "nur" "pop" "azt" "ajn" "dpl" "paj" "kla" "ptm" "bah" "btho" "hho" "pgin" "pln" "asol" "bcib" "swo" "rbl" "pci" "rrt" "bchi" "efv" "bhm" "bcm" "babt" "xop" "shq" "colw" "bmt" "pkt" "aua" "pkb" "kde" "metx" "eex" "dnx" "pmv" "fti" "eab" "kvq" "syt" "ipi" "bvm" "cpj" "lay" "pfv" "sagm" "tbd" "capn" "cko" "bama" "pol" "ptp" "ati" "ncs" "pho" "enx" "splm" "smur" "bsuv" "csh"
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