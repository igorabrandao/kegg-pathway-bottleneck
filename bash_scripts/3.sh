 
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
  "oeu" "ofo" "olu" "opf" "osa" "ota" "ots" "pace" "paj" "pam" "pamg" "pand" "paq" "pavl" "pazo" "pbar" "pbh" "pbl" "pbm" "pbn" "pca" "pcf" "pcg" "pci" "pco" "pcs" "pda" "pdag" "pdp" "pef" "pen" "peq" "pes" "peu" "pfb" "pfc" "pfe" "pff" "pfg" "pfj" "pfl" "pfn" "pfo" "pfp" "pfs" "pfv" "pfw" "pfx" "pgc" "pge" "phu" "pin" "pkt" "plf" "plg" "pman" "pmos" "pmot" "pmp" "pmu" "pmul" "pmv" "pnr" "poi" "pop" "por" "ppb" "ppd" "ppet" "ppf" "ppi" "ppj" "ppl" "ppoa" "ppp" "pprc" "ppro" "ppt" "ppuh" "ppun" "pput" "ppw" "ppx" "ppz" "prap" "psa" "psal" "psb" "psc" "psd" "psh" "psj" "psom" "psp" "pspi" "psq" "psr" "pst" "psts" "pstu" "psu" "psuw" "psy" "psyr" "psz" "pte" "ptep" "ptx" "pul" "pus" "pva" "pvm" "pvr" "pxy" "qsu" "rab" "raf" "rak" "ram" "rau" "rbd" "rbe" "rbn" "rbo" "rcc" "rcm" "rco" "rdp" "rfe" "rge" "rgl" "rhd" "rhe" "rhh" "ric" "rja" "rma" "rmc" "rmd" "rmi" "rmo" "rms" "rpg" "rph" "rpk" "rpl" "rpn" "rpo" "rpp" "rpq" "rpr" "rps" "rpv" "rpw" "rpz" "rra" "rrb" "rrc" "rre" "rrh" "rri" "rrj" "rrm" "rrn" "rrp" "rrr" "rsv" "rsw" "rtb" "rtt" "rty" "rvi" "saga" "salh" "salk" "salm" "saln" "samy" "sat" "saz" "sba" "sbb" "sbf" "sbi" "sbl" "sbm" "sbn" "sbp" "sbs" "sbt" "sbw" "scl" "scm" "scu" "sde" "sdf" "sdl" "sdn" "sdy" "sdz" "seds" "sfr" "sfu" "shal" "shd" "shm" "shp" "shq" "sind" "sita" "sku" "sla" "slh" "slim" "slo" "sly" "smo" "smul" "snn" "soc" "sod" "soe" "sok" "son" "sot" "spc" "spen" "spl" "spo" "spoi" "sros" "ssal" "sse" "stem" "stes" "sua" "sulj" "sulr" "suls" "sun" "sur" "sutk" "sutt" "sva" "tao" "tas" "tat" "tau" "tbn" "tca" "tcx" "tdn" "tea" "tee" "teg" "tgr" "thes" "thi" "thig" "thin" "thio" "tho" "tig" "tin" "tkm" "tmb" "tml" "tms" "tni" "tnl" "tsn" "tsy" "ttc" "tti" "ttu" "tve" "tvi" "tvr" "tvs" "uma" "upv" "ure" "vce" "vcf" "vch" "vci" "vcj" "vcl" "vcn" "vco" "vcq" "vcr" "vcs" "vcz" "vem" "vin" "vok" "vvi" "vvu" "wbr" "wgl" "wic" "wma" "woc" "wsu" "xba" "xbc" "ypa" "ypd" "ype" "ypg" "yph" "ypj" "ypk" "ypl" "ypm" "ypn" "ypo" "ypp" "yps" "ypt" "ypv" "ypw" "ypx" "ypz" "zal" "zdf" "zin" "zma" "zne" "ztr" "aaa" "aacn" "aact" "aaf" "aaqu" "abo" "acan" "acid" "acin" "acio" "ack" "acn" "acom" "acp" "acra" "acx" "ade" "adf" "afa" "afq" "ago" "ajs" "aln" "aly" "ama" "amah" "amf" "amil" "amim" "amp" "amw" "ank" "aoa" "aoh" "apac" "apag" "apd" "aph" "apha" "app" "apro" "aql" "aqs" "aqu" "aseg" "asy" "ath" "atw" "avd" "avl" "avn" "awd" "axe" "aza" "azi" "azo" "azr" "bbo" "bced" "bcep" "bcon" "bct" "bdf" "bdl" "beq" "bgd" "bge" "bgf" "bgl" "bgu" "blat" "bmec" "bna" "boe" "bph" "bpla" "bpsi" "bpsl" "bpx" "bpy" "bsem" "bstl" "btei" "btra" "btrh" "bub" "bud" "bue" "bug" "bui" "buk" "buo" "buq" "buz" "bxb" "cal" "camp" "caq" "caur" "cbra" "cco" "ccoc" "ccv" "ccx" "cdn" "cdu" "cfar" "cfd" "cfp" "cft" "cfx" "cfz" "cgr" "cha" "cho" "chrb" "chri" "chrm" "chro" "cic" "cif" "cir" "cit" "cjb" "cjd" "cjej" "cjen" "cjer" "cjeu" "cji" "cjj" "cjl" "cjm" "cjn" "cjp" "cjq" "cjr" "cjs" "cjv" "cjw" "cjx" "cjy" "cjz" "cko" "clu" "clus" "cmai" "cme" "cola" "colw" "com" "cot" "cov" "coz" "cpap" "cpot" "cps" "crb" "crg" "crh" "cri" "cro" "crt" "cru" "crv" "crz" "csat" "cten" "cthr" "ctp" "ctt" "cvc" "cwe" 
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