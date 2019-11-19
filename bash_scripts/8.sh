 
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
  "bmt" "sez" "bar" "lmf" "jeo" "rbg" "lpx" "bsz" "bat" "kur" "sagu" "lgy" "apf" "mtw" "lmog" "kzo" "thaa" "lfu" "apu" "pleo" "lmp" "spsy" "geh" "sezo" "lor" "lys" "bsv" "bah" "spor" "sequ" "lpar" "lyb" "bsw" "bai" "spop" "salo" "seu" "lle" "apg" "lmol" "bsg" "mey" "apq" "lmoj" "bov" "maad" "lpd" "apx" "lmoz" "bcs" "sure" "mmed" "sedi" "lyg" "bsk" "spos" "aua" "bant" "hml" "hhd" "bol" "brn" "banr" "sds" "boo" "apz" "bans" "sdg" "lcu" "pseb" "apk" "banh" "lct" "lit" "asz" "hmn" "rst" "bcar" "hci" "paek" "bcas" "hct" "hcc" "banv" "ocd" "hcd" "hli" "asv" "rbz" "mcg" "panc" "bmr" "thas" "thd" "aae" "aalg" "abi" "acia" "actt" "aev" "afd" "agi" "agl" "ahz" "ams" "amu" "amy" "ane" "ape" "aqb" "arac" "arb" "arf" "ark" "asan" "asc" "ase" "asg" "asl" "aue" "aui" "aun" "aur" "avs" "baab" "bafe" "bafh" "baft" "bana" "bane" "bbau" "bbgw" "bbs" "bchi" "bgw" "bhf" "bhp" "bih" "bip" "bis" "blb" "blj" "blo" "blp" "bmay" "bprl" "bpw" "bth" "btho" "bvs" "bvt" "caa" "caci" "cacn" "cae" "camg" "capn" "capq" "cay" "cbb" "cbf" "cbl" "cbm" "cbn" "cbol" "ccm" "cdiv" "cfi" "cfl" "cgrn" "cgw" "chz" "clg" "clp" "cmp" "coc" "col" "cpe" "cpeg" "cpf" "cpj" "cpm" "cpr" "cpro" "cprv" "cpsd" "cpst" "csh" "cso" "csph" "ctcf" "cthi" "ctrk" "dau" "dco" "dfg" "dku" "dmi" "dod" "dte" "dun" "ehl" "eol" "fac" "fae" "fat" "fcm" "fco" "fek" "ffa" "fgg" "fgl" "fgo" "fhw" "fin" "fki" "fla" "fli" "flm" "fln" "fmo" "fnc" "fne" "fnf" "fnt" "fnu" "foh" "fop" "for" "fpal" "fpd" "fpk" "fpo" "fpq" "fpw" "fsi" "fsl" "ful" "fus" "fva" "gac" "gez" "gor" "gtl" "halz" "hgi" "hho" "hhy" "hjt" "hlc" "hlr" "hme" "hru" "hsd" "hsw" "hte" "hth" "hxa" "hya" "hyd" "hye" "hys" "ipo" "kab" "kit" "kos" "ksk" "laci" "lan" "lba" "lbi" "lbj" "lbl" "lbo" "lby" "len" "leo" "leq" "let" "lhf" "lic" "lie" "lil" "lis" "lot" "lpil" "lua" "lul" "lut" "lvn" "mab" "marc" "marf" "max" "mbg" "mbu" "mbw" "mby" "mcb" "mcj" "mcn" "mcw" "mear" "meg" "mel" "mem" "mema" "mer" "mev" "mfz" "mgan" "mge" "mgf" "mgh" "mgi" "mgin" "mgn" "mgt" "mgv" "mgw" "mhaz" "mhi" "mhor" "mhu" "mhz" "mib" "micb" "mich" "micr" "mil" "min" "mip" "mja" "mkc" "mla" "mli" "mll" "mlt" "mmae" "mmet" "mmh" "mmi" "mod" "mpal" "mpi" "mpru" "mpy" "mrh" "mro" "mrs" "mse" "mthe" "mthn" "mthr" "mtp" "muh" "mvq" "myv" "mzh" "nac" "nag" "naj" "nan" "nat" "nbg" "ndo" "nfa" "nge" "nin" "nmg" "nob" "noj" "nom" "nou" "npe" "npl" "nsd" "nvr" "oaq" "ock" "oll" "orh" "ori" "pab" "pac" "pach" "pacn" "pah" "parc" "pary" "pbor" "pbt" "pcre" "pdi" "pep" "pet" "pfi" "pfr" "pfre" "pfu" "pgi" "pgin" "pgn" "pgo" "pgt" "phe" "pho" "pko" "pmn" "pmuc" "pmx" "pob" "pom" "ppn" "prf" "prl" "pseg" "pseu" "psez" "pth" "pto" "ptq" "puf" "pwo" "pya" "pyc" "pyn" "pys" "rae" "rhu" "ria" "roa" "rop" "rpla" "rpy" "rqi" "rrt" "rrz" "rsd" "rtc" "rtn" "run" "rup" "saal" "saci" "saf" "sali" "salt" "sbag" "scj" "sck" "scla" "scou" "sele" "selo" "selt" "seon" "sfk" "sfz" "sge" "sgn" "shj" "ske" "skn" "sll" "slp" "slr" "slw" "slx" "smr" "snz" "soo" "sox" "spit" "spon" "sprn" "srb" "srg" "sri" "ssg" "sth" "stro" "sul" "swo" "syc" "syf" "syi" "syj" "syl" "synw" "syo" "syq" "sys" "syt" "syw" "tac" "taj" "tal" "tam" "tar" "tbw" "tcb" "tdi" "ten" "tfo" "thb" "thf" "tje" "tko" "tma" "tmar" "tmg" "tmi" "tmm" "tmq" "toh" "tpe" "tprf" "trd" "ttk" "tvo" "udi" "vai" "vba" "vbs" "vdi" "vma" "vmo" "wba" "wcf" "wfu" "wpa" "wso" "wwe" "xii" "xyl" "zga"
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