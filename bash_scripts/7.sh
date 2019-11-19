 
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
  "spo" "sdt" "mvq" "run" "sphy" "ecoa" "mza" "hsm" "dso" "bko" "rhp" "afr" "abq" "meta" "sag" "sya" "hmg" "bme" "pig" "sod" "mkc" "spn" "hpaz" "cjeu" "babs" "ftq" "ham" "kmx" "tml" "mah" "csa" "eba" "nge" "pfu" "spen" "gxy" "rua" "stv" "ppla" "apg" "lus" "ffa" "kpk" "apx" "sauu" "psj" "seg" "cek" "bge" "pst" "lie" "atd" "spv" "ypt" "woc" "bpt" "bdf" "cae" "cjx" "bcd" "lac" "bna" "ctp" "ptu" "mgan" "dvg" "bmae" "ypp" "aqs" "eun" "bon" "dds" "plq" "rab" "wcf" "cpw" "tnl" "hpyb" "idt" "erc" "rle" "myv" "bbay" "avd" "dku" "camg" "rhi" "serf" "cwe" "azi" "exf" "fli" "ypv" "mgt" "pvr" "spsw" "mli" "prap" "mlt" "api" "tin" "bpd" "salh" "seb" "blas" "cal" "plen" "stn" "ppl" "edh" "bmel" "hpyi" "aea" "tht" "pmn" "paur" "deo" "sply" "sud" "bli" "bng" "pseu" "cbl" "ptt" "pfi" "pdr" "pen" "pcz" "rlu" "pld" "pco" "qsu" "lot" "pgb" "cjn" "pge" "psw" "acio" "gez" "mcys" "spmi" "hch" "sscz" "zpa" "ask" "rrj" "fts" "shv" "ebr" "baus" "mzh" "lhv" "syj" "ecle" "sauf" "lce" "rpg" "rrh" "sagq" "lcb" "aql" "dye" "dji" "pamn" "elo" "ttc" "coc" "bay" "mthd" "mcak" "bpsa" "dea" "trd" "bgj" "wcl" "saub" "mets" "elf" "ypi" "msu" "sey" "bpsl" "acn" "rep" "fop" "bmet" "cell" "dfn" "obo" "bsup" "smr" "boh" "saly" "malg" "rak" "apu" "lmf" "bmb" "sha" "vcr" "pfc" "eclg" "psk" "yee" "chw" "leq" "les" "bmf" "pow" "eoc" "mve" "lku" "brn" "aeu" "dtr" "meg" "buf" "laqu" "spyh" "tps" "lrz" "pgz" "san" "cbd" "pay" "shd" "lpd" "ksk" "npa" "slu" "mros" "kpj" "xcw" "cmik" "vgo" "pant" "eok" "asan" "pfl" "eat" "chj" "vpf" "gbm" "spk" "bud" "cbs" "hpyu" "bacw" "cdiv" "lar" "spse" "sauz" "lpz" "rsv" "cjej" "hiw" "lpap" "cjj" "pstg" "pck" "snd" "amu" "gai" "ago" "sds" "kro" "psen" "smd" "kyr" "bdl" "rbh" "vcs" "fgg" "tsy" "lhh" "pmu" "xoy" "dno" "ssil" "das" "ctt" "cji" "btr" "ppeo" "nmc" "mxa" "noj" "stra" "ddc" "marc" "sig" "mbl" "lpl" "hmr" "amb" "hna" "ent" "thb" "cjen" "spkc" "vta" "ypw" "mai" "het" "elg" "bov" "bmy" "net" "hpak" "atf" "xcb" "neo" "rbu" "rhr" "clx" "noc" "bor" "wol" "toh" "hsw" "phc" "sak" "adt" "rfe" "maad" "rtt" "fco" "lgn" "bafh" "clp" "fsl" "bcj" "stg" "caci" "mhq" "sch" "mtim" "ahn" "mtw" "clus" "ahz" "cho" "ndo" "ema" "bww" "mfz" "htr" "clec" "spl" "bcc" "aaj" "mya" "cft" "got" "sepp" "sva" "hic" "gox" "babo" "cohn" "spq" "ppd" "ekf" "pse" "kpg" "lpq" "acx" "psos" "asz" "seu" "fnu" "mhu" "ane" "banv" "aui" "blh" "eic" "mrh" "aly" "lmp" "ptc" "sea" "soc" "ebh" "hli" "apw" "sagp" "muh" "sulj" "pput" "pch" "jcu" "pste" "bsh" "afe" "rhd" "nou" "fiy" "rpn" "kbl" "oar" "lpj" "nhl" "cpro" "fth" "sri" "bur" "lmu" "tfo" "nba" "xcp" "pgt" "xcu" "csat" "ddi" "ggh" "tho" "sty" "ztr" "bhp" "ehd" "tve" "baa" "bapi" "kpe" "bit" "hhs" "ssuy" "cey" "mmed" "pmuc" "aaw" "bmz" "ssb" "rvi" "pdp" "bhf" "sdg" "bmaq" "amyt" "xca" "zga" "afj" "ncr" "nmt" "bmi" "sulf" "hth" "bamb" "bau" "vex" "eclz" "pavl" "dfa" "ppt" "ppol" "spj" "otm" "cpap" "pmp" "gya" "wsu" "msd" "sauw" "gsu" "tmz" "slp" "panc" "hsz" "htl" "adl" "thas" "abf" "chrb" "miv" "rsl" "mgl" "cbf" "ppo" "cdq" "jeo" "miu" "tci" "psb" "hso" "ban" "ipo" "ctu" "aseg" "sjp" "tcl" "dgg" "tmq" "agn" "bpl" "hpd" "lrl" "suli" "cgi" "ena" "gdj" "sxl" "lvn" "buq" "bper" "ptl" "mrm" "pcx" "bsq" "ypc" "oco" "ori" "mtp" "pfz" "rrb" "xcn" "lmt" "dpp" "rhf" "rge" "ddn" "chz" "fnc" "bby" "paeq" "ypk" "bst" "our" "mhz" "pseb" "bms" "sage" "pkr" "mep" "smy" "pgr" "ljh" "bmal" "ebi" "btab" "amw" "gtm" "shn" "sfz" "mhao" "mrr" "lez" "enm" "mgi" "mpr" "pva" "rco" "zin" "dpg" "esi" "mbr" "apj" "laz" "agi" "azr" "wbm" "nve" "aht" "rsm" "mmet" "hbl" "bip" "pace" "psil" "erj" "ppn" "cfu" "gps" "lax" "fat" "lei" "vpd" "rlg" "vma" "tvs" "vsc" "aae" "sps" "banc" "pfg" "ppq" "gan" "echv" "bbr" "gho" "egn" "meu" "rhc" "tad" "hpye" "bhz" "suh" "ssy" "sxo" "stz" "aci" "aza" "lmir" "lhk" "pan" "bav" "bsy" "bpdz" "bsz" "pseo" "ocm" "lmo" "man" "ssg" "aamy" "ttt" "spha" "beo" "lyb" "asl" "tpty" "masw" "bpx" "haf" "baro" "ypb" "cel" "pko" "jab" "tma" "mhx" "baj" "acep" "gjf" "elw" "doa" "sbn" "mcl" "jeu" "pac" "hay" "azz" "vsp" "sthe" "bapg" "kak" "lhu" "rez" "oaq" "tpq" "epe" "sug" "tcb" "gbe" "sphu" "sab" "echl" "don" "ots" "hey" "yrb" "vvl" "mfa" "esl" "nat" "mthn" "hmn" "xbo" "kpt" "len" "wba" "aep" "scon" "tgo" "ccon" "bat" "stu" "rps" "srq" "paf" "meti" "lrr" "agr" "aat" "fus" "bcom" "bop" "koe" "hem" "pha" "pic" "aqb" "tac" "eli" "hen" "ssf" "mmas" "bcep" "bgf" "bbar" "bgw" "buz" "aeh" "vaa" "bpe" "kpw" "ppy" "pmul" "saui" "paen" "nri" "cmar" "ftv" "banr" "lan" "dsa" "arv" "hcp" "sauv" "hhm" "crb" "hml" "gbs" "panr" "syi" "thw" "prl" "fin" "csph" "srg" "tkm" "ykr" "ebd" "ncb" "six" "ccav" "izh" "ppb" "lcd" "rso" "pswu" "bck" "acv" "ndi" "mec" "mbah" "spit" "bgy" "axo" "xbc" "rtc" "pfre" "ysi" "sulz" "ljf" "blag" "lrm" 
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