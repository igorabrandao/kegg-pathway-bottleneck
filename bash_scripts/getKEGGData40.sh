 
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
  "bde" "bdn" "bbi" "bbp" "bbf" "bbv" "bbru" "bbre" "bbrv" "bbrj" "bbrc" "bbrn" "bbrs" "bbrd" "bast" "btp" "bcor" "bka" "bks" "bcat" "bpsp" "bii" "bang" "bpsc" "bsca" "bact" "bcho" "gva" "gvg" "gvh" "sij" "pdo" "tbi" "aey" "plan" "plak" "plim" "psuf" "pvs" "pvn" "abam" "nhi" "nab" "abai" "rxy" "rrd" "bsol" "cwo" "afo" "aym" "ccu" "shi" "ele" "eyy" "gpa" "aeq" "ddt" "apv" "ols" "olo" "pcat" "cgo" "caer" "cbac" "erz" "euz" "syn" "syz" "syy" "syt" "sys" "syq" "syj" "syo" "syw" "syc" "syf" "syd" "sye" "syg" "syr" "syx" "syp" "cya" "cyb" "syne" "synp" "synk" "synr" "synd" "syu" "syh" "synw" "slw" "syv" "syl" "tel" "thn" "tvn" "thec" "cgc" "cyi" "dsl" "cmp" "lep" "len" "let" "lbo" "hhg" "pseu" "pser" "pma" "pmm" "pmt" "pmn" "pmi" "pmb" "pmc" "pmf" "pmg" "pmh" "pmj" "pme" "prc" "prm" "amr" "glp" "gen" "gee" "chon" "mar" "mpk" "miq" "mvz" "can" "csn" "cyl" "hao" "enn" "cyu" "cyt" "cwa" "cyp" "cyh" "cyc" "cyj" "cyn" "ter" "mic" "arp" "pagh" "gei" "oac" "oni" "mpro" "cep" "gvi" "glj" "ana" "npu" "nos" "nop" "non" "nfl" "noe" "ava" "naz" "anb" "acy" "awa" "csg" "calo" "calt" "calh" "riv" "fis" "nsp" "dou" "ncn" "cthe" "plp" "scs" "stan" "ceo" "cer" "mbf" "det" "deh" "deb" "dev" "deg" "dmc" "dmd" "dmg" "dmx" "dmy" "dmz" "duc" "dly" "dew" "dfo" "rrs" "rca" "cau" "chl" "cag" "hau" "tro" "sti" "atm" "abat" "psub" "abao" "cap" "pbf" "kbs" "dra" "dge" "ddr" "dmr" "dpt" "dgo" "dpd" "dsw" "dch" "dab" "dpu" "dez" "dwu" "dfc" "dein" "tra" "tth" "ttj" "tts" "ttl" "tsc" "thc" "tos" "taq" "tpar" "tbc" "mrb" "mre" "msv" "mtai" "opr" "mhd" "ccz" "fgi" "ttr" "ctr" "ctd" "ctf" "ctrd" "ctro" "ctrt" "cta" "cty" "cra" "ctrq" "ctrx" "ctrz" "ctrp" "ctlj" "ctlx" "ctll" "ctb" "ctrr" "ctlf" "ctli" "ctl" "ctru" "ctrl" "ctrv" "ctrm" "ctla" "ctlm" "ctls" "ctlz" "ctlc" "ctln" "ctlb" "ctlq" "cto" "ctrn" "ctj" "ctz" "ctg" "ctk" "csw" "ces" "ctrb" "ctre" "ctrs" "ctec" "cfs" "cfw" "ctfw" "ctrf" "ctch" "ctn" "ctq" "ctv" "ctw" "ctrg" "ctri" "ctra" "ctrh" "ctrj" "ctrk" "ctjt" "ctcf" "ctfs" "cthf" "ctcj" "cthj" "ctmj" "cttj" "ctjs" "ctrc" "ctrw" "ctry" "ctct" "cmu" "cmur" "cmn" "cmm" "cmg" "cmx" "cmz" "cpn" "cpa" "cpj" "cpt" "clp" "cpm" "cpec" "cpeo" "cper" "chp" "chb" "chs" "chi" "cht" "chc" "chr" "cpsc" "cpsn" "cpsb" "cpsg" "cpsm" "cpsi" "cpsv" "cpsw" "cpst" "cpsd" "cpsa" "cav" "cca" "cab" "cabo" "cfe" "cgz" "chla" "pcu" "pnl" "puv" "ney" "wch" "sng" "ote" "obg" "vbh" "obt" "caa" "amu" "agl" "xii" "min" "mkc" "vba" "vbs" "rba" "psl" "pir" "plm" "peh" "pbs" "pls" "plh" "fmr" "ttf" "rul" "gmr" "mff" "ges" "gog" "ipa" "saci" "pbor" "agv" "kst" "phm" "pbu" "pbp" "pbas" "vbl" "vbc" "bbu" "bbz" "bbn" "bbj" "bbur" "bga" "bgb" "bgn" "bgs" "bgc" "baf" "bafz" "bafh" "baft" "bafe" "bbs" "bvt" "bchi" "bmay" "btu" "bhr" "bhi" "bdu" "bre" "bcw" "bmo" "bmiy" "bpak" "bane" "btur" "tpa" "tpw" "tpp" "tpu" "tph" "tpo" "tpas" "tpc" "tpg" "tpm" "tpb" "tde" "tsu" "tbe" "taz" "tpi" "tpl" "tped" "scd" "tpk" "trm" "ssm" "sta" "stq" "sfc" "sper" "sbu" "scc" "sgp" "slr" "ock" "lil" "lie" "lic" "lis" "lbj" "lbl" "lbi" "lbf" "lst" "laj" "lmay" "lkm" "lwl" "tpx" "bhy" "bhd" "brm" "bpo" "bpj" "bpip" "bpw" "bip" "bhp" "aba" "aca" "acm" "gma" "tsa" "trs" "talb" "abas" "sus" "ctm" "abac" "emi" "epo" "eti" "rsd" "fnu" "fnc" "fnt" "fus" "fne" "fhw" "fpd" "fva" "ful" "fmo" "fgo" "fnf" "ipo" "lba" "leo" "lot" "leq" "lhf" "str" "smf" "sns" "tai" "aco" "tli" "amo" "sbr" "cpor" "fsu" "fsc" "gau" "gph" "gba" "bth" "btho" "bfr" "bfs" "bfg" "bfb" "bvu" "bhl" "bsa" "bxy" "bdo" "bdh" "boa" "bcel" "bcac" "bcae" "bzg" "bhf" "bis" "pgi" "pgn" "pgt" "pah" "pcre" "pbt" "pmuc" "pet" "ppn" "pdi" "parc" "tfo" "toh" "pary" "dun" "bvs" "psac" "osp" "buy" "aps" "pru" "pmz" "pdn" "pit" "pdt" "pro" "pfus" "peo" "pje" "poc" "alq" "afd" "ash" "rbc" "bacc" "dori" "blq" "asx" "mbas" "sru" "srm" "rmr" "rmg" "rbar" "cpi" "cbae" "chit" "nko" "nso" "nia" "fla" "fgg" "arb" "ark" "agi" "arac" "fln" "pseg" "pgin" "pgo" "hhy" "sgn" "phe" "pep" "pcm" "psty" "pgs" "psn" "shg" "sht" "sphn" "smiz" "spsc" "sphz" "scn" "mup" "muc" "mgot" "muh" "mgin" "mgk" "agd" "oli" "sbx" "cmr" "camu" "bbd" "evi" "est" "echi" "alm" "chu" "dfe" "sli" "srd" "smon" "spir" "spik" "lby" "rsi" "run" "rup" "eol" "fae" "fib" "psez" "als" "fli" "hsw" "hym" "hyd" "hye" "hyg" "hyp" "hyz" "hnv" "hyh" "hyj" "pko" "pact" "ruf" "rti" "rud" "mtt" "flm" "fll" "fpf" "fbt" "aas" "che" "cec" "cher" "chk" "gfo" "grl" "gfl" "grs" "fjo" "fjg" "fps" "fpc" "fpy" "fpo" "fpq" "fpv" "fpw" "fpk" "fpsz" "fbr" "fco" "fin" "fgl" "fcm" "ffa" "fat" "fki" "fpal" "coc" "ccm" "col" "chg" "capn" "cgh" "clk" "cspu" "ccyn" "caph" "csto" "capq" "rbi" "zpr" "cat" "ran" "rai" "rar" "rag" "rae" "rat" "fbc" "marm" "mart" "marb" "mare" "cao" "cly" "clh" "cbal" "cbat" "wvi" "kdi" "dok" "ddo" "dod" "lan" "lvn" "laci" "zga" "mrs" "mlt" "asl" "aev" "orh" "ori" "ptq" "ndo" "nom" "nsd" "nob" "noj" "pom" "pob" "prn" "pola" "poa" "eao" "emn" "een" "elb" "emg" "ego" "egm" "elz" "myr" "mpw" "mod" "myz" "chz" "cgn" "cih" "chh" "cio" "chry" "cpip" "ctak" "chrs" "chrz" "carh" "csha" "win" "wij" "sze" "ahz" "syi" "tdi" "ten" "tje" "tmar" "lut" "lul" "wfu" "for" "foh" "fop" "salt" "seon" "aalg" "oll" "oaq" "fek" "taj" "aue" "spon" "kos" "marf" "aqb" "aqa" "aqd" "emar" "cnr" "mur" "psyn" "afla" "anp" "ebv" "fba" "fbu" "fbe" "smg" "sms" "smh" "sum" "smv" "smub" "smum" "smue" "smup" "bbl" "bpi" "bmm" "bcp" "bbg" "bbq" "blp" "blu" "blck" "fte" "flu" "oho" "ise" "elv" "udi" "bbau" "cte" "cpc" "clz" "cch" "cph" "cpb" "cli" "pvi" "plt" "pph" "paa" "proc" "prs" "pros" "cts" "ial" "mro" "cprv" "caci" "aae" "hya" "hho" "hys" "hth" "hte" "tal" "trd" "sul" "saf" "pmx" "ttk" "tam" "dte" "tma" "tmm" "tmi" "tmw" "tmq" "tmx" "tpt" "trq" "tna" "tnp" "thq" "thz" "thr" "tle" "tta" "phy" "tme" "taf" "thp" "ther" "fno" "fpe" "fia" "pmo" "mpz" "marn" "dtn" "kol" "kpf" "mpg" "minf" "cex" "din" "ddf" "dap" "cni" "fsi" "gtl" "caby" "dth" "dtu" "tye" "nde" "nmv" "nio" "nja" "lfc" "lfi" "lfp" "leg" "tid" "top" "tcm" "thet" "cthi" "saal" "sbe" "sbag" "sox" "prf" "bana" "bih" "vai" "mox" "dpb" "tmg" "srb" "srg" "wwe" "bgw" "bbgw" "mib" "wba" "pwo" "cgw" "baab" "mja" "mfe" "mvu" "mfs" "mif" "mjh" "mig" "mmp" "mmq" "mmx" "mmz" "mmd" "mmak" "mmao" "mae" "mvn" "mvo" "mok" "metf" "mth" "mmg" "metc" "mwo" "mete" "metz" "mst" "metb" "msi" "mru" "meb" "mmil" "meye" "mol" "mel" "mew" "meth" "mfc" "mfi" "mcub" "msub" "metn" "mett" "meto" "mfv" "mka" "afu" "afg" "apo" "ave" "ast" "fpl" "gac" "gah" "tac" "tvo" "pto" "fac" "fai" "cdiv" "tar" "max" "mer" "mear" "marc" "abi" "acf" "pho" "pab" "pfu" "pfi" "pyn" "pya" "pys" "pyc" "tko" "ton" "tga" "tsi" "tba" "the" "tha" "thm" "tlt" "ths" "tnu" "teu" "tgy" "thv" "tch" "tpep" "tpie" "tgg" "tce" "tbs" "thh" "tsl" "ttd" "tprf" "trl" "tpaf" "thy" "ppac" "mac" "mba" "mby" "mbw" "mbar" "mbak" "mma" "mmaz" "mmj" "mmac" "mvc" "mek" "mls" "metm" "mef" "meq" "msj" "msz" "msw" "mthe" "mthr" "mhor" "mfz" "mbu" "mmet" "mmh" "mhaz" "mev" "mzh" "mpy" "mhz" "mtp" "mcj" "mhi" "mhu" "mla" "mem" "mbg" "mema" "mpi" "mbn" "mfo" "mpl" "mpd" "mez" "rci" "hal" "hsl" "hdl" "hhb" "hje" "halh" "hhsr" "hsu" "hsf" "salr" "hma" "hhi" "hhn" "hab" "hta" "nph" "nmo" "hut" "hti" "hmu" "halz" "hali" "hsn" "harc" "hwa" "hwc" "hvo" "hme" "hgi" "hbo" "haq" "haj" "haer" "hlm" "halm" "hla" "halp" "halb" "hezz" "srub" "hae" "hah" "htu" "hda" "hjt" "nmg" "hxa" "nat" "npe" "nvr" "npl" "nge" "hru" "nou" "sali" "hlr" "hlc" "naj" "nag" "nan" "nbg" "nac" "ape" "acj" "smr" "shc" "iho" "iis" "dka" "dfd" "dmu" "tag" "iag" "thg" "hbu" "pfm" "pdl" "sto" "sso" "sol" "ssoa" "ssol" "ssof" "sai" "sacn" "sacr" "sacs" "sis" "sia" "sim" "sid" "siy" "sin" "sii" "sih" "sir" "sic" "sula" "mse" "mcn" "mhk" "mpru" "aho" "aman" "abri" "asul" "sacd" "pai" "pis" "pcl" "pas" "pyr" "pog" "tne" "pyw" "cma" "tuz" "ttn" "vdi" "vmo" "tpe" "thb" "tcb" "thf" "asc" "acia" "clg" "ffo" "nmr" "nir" "nkr" "nid" "nin" "niw" "nct" "csy" "nga" "nvn" "nev" "taa" "nfn" "ncv" "csu" "nbv" "tah" "ndv" "neq" "naa" "marh" "kcr" "barc" "barb" "loki" "agw" "arg"
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