 
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
  "abaa" "abau" "achb" "adt" "afm" "ahh" "ahi" "ahj" "aka" "amac" "amar" "amb" "ani" "aon" "asw" "axn" "axo" "axy" "bab" "bac" "baj" "bajc" "bak" "bap" "bapf" "bapg" "baph" "bapu" "bapw" "bas" "bau" "bav" "baw" "bba" "bbac" "bbag" "bbat" "bbay" "bbh" "bbm" "bbr" "bbro" "bbw" "bcc" "bcen" "bceo" "bcew" "bch" "bci" "bcib" "bcig" "bcj" "bcm" "bcn" "bcom" "bdi" "bdq" "bex" "bfz" "bgj" "bgm" "bhm" "bho" "bhz" "blep" "bma" "bmab" "bmae" "bmaf" "bmai" "bmal" "bman" "bmaq" "bmaz" "bmj" "bmk" "bml" "bmn" "bmor" "bmu" "bmul" "bmv" "bmx" "bmy" "bng" "boc" "boh" "bok" "bor" "boz" "bpa" "bpc" "bpd" "bpdz" "bpe" "bper" "bpet" "bpeu" "bpg" "bpk" "bpl" "bpm" "bpq" "bpsa" "bpsd" "bpse" "bpsh" "bpsm" "bpso" "bpsu" "bpt" "bpz" "brb" "bsc" "bsed" "bsto" "btab" "btd" "bte" "bter" "btha" "bthe" "bthl" "bthm" "btj" "btq" "btrm" "btv" "btz" "bua" "buc" "buf" "bup" "bur" "but" "buu" "bve" "bvg" "bvi" "bze" "camy" "camz" "cann" "care" "cate" "cavi" "cbc" "cbd" "cbg" "cbr" "cbs" "cbu" "ccav" "cci" "ccon" "ccun" "cdm" "cea" "ced" "cek" "cel" "cell" "cey" "cfo" "cfu" "cfus" "cgi" "chj" "chw" "cim" "cja" "clec" "clx" "cmed" "cmik" "cmj" "cmw" "cne" "cox" "cpaf" "cput" "cpw" "cqi" "csa" "csi" "csj" "csk" "csl" "csol" "cspf" "csz" "ctu" "cui" "cvr" "cyq" "cyy" "dao" "dav" "dax" "dcr" "dct" "dda" "ddc" "ddn" "ddq" "dds" "den" "des" "deu" "dfn" "dic" "dja" "dji" "dko" "dno" "dnx" "dosa" "dpa" "dpl" "dpte" "dpx" "dso" "dsq" "dti" "dtx" "dvg" "dvl" "dvu" "dye" "dzc" "dze" "eab" "eam" "eas" "eau" "eay" "ebc" "ebd" "ebe" "ebf" "ebh" "ebi" "ebl" "ebr" "ebs" "ebu" "ecan" "ecg" "eci" "eck" "ecl" "ecla" "eclc" "ecle" "eclg" "eclo" "ecls" "ecly" "eclz" "ecm" "ecoa" "ecoh" "ecol" "ecp" "ecq" "ecr" "ect" "ecv" "ecw" "ecx" "ecy" "ecz" "edh" "edj" "eec" "ege" "egu" "ehd" "ehm" "eih" "ekb" "ekf" "eko" "elc" "eld" "elf" "elg" "elh" "ell" "eln" "elo" "elp" "elr" "elu" "elw" "ena" "enc" "enl" "enm" "eno" "enr" "ent" "enx" "eoc" "eof" "eoh" "eok" "epe" "epr" "erj" "ern" "esa" "esc" "ese" "esl" "eso" "eta" "eum" "eun" "exf" "fad" "fau" "fbl" "fcd" "fgu" "fmi" "fpp" "frc" "frm" "gag" "gai" "gao" "gap" "gbi" "gbm" "geb" "ged" "gem" "glo" "gme" "gni" "gpb" "gpi" "gps" "gsb" "gsk" "gsu" "gtr" "gur" "hac" "had" "hahe" "hak" "halo" "ham" "han" "hap" "hax" "hay" "haz" "hbi" "hbl" "hbr" "hcb" "hce" "hch" "hcn" "hco" "hcp" "hcs" "hdu" "heb" "hee" "hef" "hei" "hel" "hem" "hen" "hep" "hes" "het" "heu" "hex" "hey" "hfe" "hha" "hhc" "hhe" "hhf" "hhk" "hhm" "hht" "hhu" "hhz" "hia" "hic" "hih" "hik" "hir" "hiw" "hix" "hja" "hmi" "hmr" "hms" "hna" "hoh" "hpak" "hpas" "hpaz" "hpd" "hpf" "hph" "hpit" "hpn" "hpr" "hpt" "hpx" "hpya" "hpyb" "hpyc" "hpyd" "hpye" "hpyf" "hpyg" "hpyh" "hpyi" "hpyj" "hpyk" "hpym" "hpyo" "hpyr" "hpys" "hpyu" "hpz" "hrb" "hse" "hsi" "hsm" "hso" "hst" "hsz" "htr" "hty" "hyf" "hyo" "icp" "ini" "isc" "izh" "jab" "jaj" "jal" "jaz" "jcu" "jre" "kgo" "kie" "kle" "lab" "lal" "laq" "lax" "laz" "lbc" "lch" "lcp" "lee" "lef" "leh" "lei" "lew" "lez" "lgu" "lhu" "llu" "lmir" "lni" "lph" "lpn" "lpo" "lpu" "lpv" "lpy" "lrz" "lsv" "lue" "lug" "lum" "lus" "lyj" "maga" "mali" "mass" "masw" "masy" "masz" "mbah" "mbd" "mbl" "mcat" "mcs" "mct" "mcys" "mdl" "mec" "mej" "mesc" "mety" "mhae" "mhal" "mham" "mhao" "mhaq" "mhat" "mhay" "mhq" "mht" "mhx" "micc" "mii" "mis" "miu" "mmaa" "mng" "mnr" "moi" "mpha" "mpp" "mpr" "mpur" "mrm" "mrr" "msu" "mthd" "mtim" "mum" "mus" "mve" "mvg" "mvi" "mvr" "mvs" "mya" "nam" "nap" "nau" "nco" "net" "neu" "nfi" "nhl" "nii" "nis" "nit" "nlc" "nme" "nmh" "nmm" "nmp" "nms" "nmu" "noc" "nok" "npa" "nsa" "nsy" "nta" "nto" "ntt" "nur" "nvi" "nvl" "nwa" "nwr" "obo" "obr" "ocm" 
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