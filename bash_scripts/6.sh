 
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
  "cim" "cjl" "csk" "cbra" "bts" "lwi" "mcg" "tmg" "paeh" "lpc" "frc" "hya" "azc" "parc" "naj" "csz" "bma" "cpr" "mear" "asv" "nom" "bapw" "plx" "gwc" "lgl" "ell" "micb" "eno" "sali" "def" "rue" "ric" "aay" "lro" "apq" "bui" "ssw" "bpu" "pyn" "gor" "ppi" "lhl" "mpo" "rty" "ret" "lmol" "fpk" "ssan" "geq" "bsv" "vnl" "bamp" "ssi" "bjs" "acia" "pta" "smc" "ese" "bpum" "puf" "maga" "paru" "fpw" "dfl" "ged" "lbo" "ppf" "slr" "pfx" "smac" "psr" "alr" "snb" "pbj" "hex" "mcw" "nse" "skr" "sscu" "sinb" "amy" "cten" "lmoz" "syq" "pea" "pvl" "saux" "ccm" "hxa" "alv" "shu" "gtl" "bmur" "por" "pau" "sry" "kpne" "kga" "bamf" "tvi" "dtx" "cps" "anx" "cbb" "bmay" "avo" "boe" "lbj" "smav" "suo" "abaa" "jal" "aon" "myi" "swd" "parr" "hac" "sph" "lhr" "ssui" "pam" "ilo" "lys" "dic" "bhk" "mus" "stc" "nlc" "ecly" "hye" "gsb" "prt" "bpsu" "bbgw" "lyj" "cir" "pom" "smeg" "lgi" "but" "lae" "cdn" "apd" "ful" "avn" "axl" "max" "btj" "cmai" "bis" "rpp" "sent" "rph" "nsd" "sclo" "cjz" "red" "egr" "lbi" "ecg" "bth" "lhd" "pxy" "ecm" "yet" "psal" "psv" "bte" "bup" "gsr" "cbol" "cjv" "vin" "shyd" "nms" "bsx" "aal" "amil" "lat" "sba" "lle" "lpo" "lib" "sok" "cro" "slh" "com" "xfh" "sup" "lel" "pdam" "vce" "avs" "ams" "lmoe" "agc" "tdn" "crz" "buo" "gag" "bpsm" "cyy" "adf" "bua" "pya" "rtr" "rhz" "ags" "xcj" "pprf" "shm" "cdu" "camp" "mhor" "crh" "bhe" "rli" "nvi" "flm" "spb" "abv" "bacs" "aro" "pli" "pia" "vha" "spc" "eci" "jaj" "pkz" "vcz" "pth" "bwe" "son" "ofo" "bman" "chri" "apac" "sgn" "slx" "hte" "gmx" "pcr" "apz" "esa" "boi" "pbr" "vba" "tmb" "siq" "pach" "hdu" "bthl" "spx" "ajs" "xnm" "sys" "ery" "pkn" "sfd" "ank" "bcn" "hbr" "cann" "ass" "had" "spat" "sams" "send" "bql" "err" "eof" "ppro" "bsto" "rir" "paek" "eclo" "pes" "psa" "eas" "fae" "spp" "bsaf" "bacb" "seon" "lhe" "xpo" "fac" "mdl" "fnt" "ppuu" "ecx" "nwi" "ttu" "camy" "psez" "rja" "sew" "sca" "hpx" "shp" "paly" "lpv" "pnr" "strn" "naf" "bbev" "shal" "cja" "fpq" "psz" "boc" "gmo" "mis" "pmw" "pans" "thu" "bmj" "bacq" "tko" "pda" "bag" "shh" "hpt" "nap" "marf" "ypa" "aph" "acd" "pfj" "pprc" "panp" "mhay" "rhk" "ala" "gse" "yhi" "abi" "saal" "awd" "babc" "rhv" "lmoc" "hpit" "salk" "nel" "nek" "ctcf" "slo" "chro" "acum" "pih" "snm" "eme" "maes" "lug" "afv" "pgd" "metr" "gmc" "olu" "esc" "rcm" "den" "upv" "snn" "salo" "yew" "lmom" "lga" "nsa" "plg" "bze" "sags" "mgin" "bsp" "des" "sagn" "npe" "sbm" "mhae" "acc" "phu" "fmo" "xyl" "bane" "hph" "bmyo" "bced" "snp" "fau" "lch" "stm" "hhk" "vam" "sutt" "nte" "ecoh" "pcre" "arf" "pmos" "mmi" "cne" "snl" "stro" "pspg" "pazo" "gac" "setu" "srp" "mpy" "awu" "zne" "psy" "suy" "slw" "acan" "csi" "ure" "zma" "bpg" "plf" "ste" "ecoj" "msc" "cci" "seen" "cmw" "las" "suls" "mby" "nmg" "bwh" "heu" "fgu" "aact" "wfu" "stes" "spne" "phe" "afm" "cij" "rst" "mbur" "epr" "apf" "sei" "rdp" "apor" "lpp" "blep" "lut" "bbw" "ang" "rsh" "dcr" "btd" "daa" "rsz" "mbu" "dsf" "asr" "anh" "ptn" "cpeg" "tvo" "hlr" "sutk" "prw" "syb" "cfp" "ppp" "pnp" "psyg" "vfi" "metl" "pod" "fpo" "xba" "bqr" "sers" "syo" "ypo" "cola" "sdf" "mno" "mgf" "sagu" "pbb" "slj" "xcm" "bbau" "hys" "kpnu" "sks" "ahi" "scm" "ssp" "nii" "thf" "pvg" "syf" "kmi" "law" "vff" "ekb" "pib" "bbac" "aka" "amf" "hmi" "wbr" "lagg" "gej" "mich" "rau" "saut" "syw" "mor" "soo" "bvt" "mmai" "psed" "bmaf" "cpot" "pmut" "barj" "dat" "abau" "rre" "dps" "xct" "slq" "hhd" "capq" "lko" "smm" "seni" "crg" "fto" "sip" "apag" "mca" "ale" "mpp" "mll" "vau" "bsg" "cthi" "macr" "cyo" "pleo" "apb" "bmyc" "caur" "rpz" "cfz" "weo" "chrm" "lhf" "sagc" "sagr" "mmn" "spon" "rpl" "hpyd" "bthm" "lpar" "pfo" "for" "php" "deu" "sulr" "tdl" "nhm" "elh" "bsm" "hmc" "sphi" "amg" "tpe" "cco" "atu" "ade" "dti" "cna" "lpy" "bvg" "banh" "hlc" "roa" "phq" "saud" "sagg" "ehi" "thaa" "sux" "hsi" 
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