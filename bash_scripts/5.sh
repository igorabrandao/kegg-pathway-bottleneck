 
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
  "ptm" "ptn" "ptt" "ptu" "pur" "pvg" "pvl" "pvu" "pvx" "pyo" "rac" "rbh" "rbu" "rfr" "rhf" "rsc" "rse" "rsl" "rsm" "rso" "rsz" "salv" "saly" "sbj" "sdr" "sea" "seb" "sec" "see" "seeh" "seen" "sef" "seg" "sega" "seh" "sei" "sek" "sel" "sena" "send" "seni" "senl" "senn" "seno" "senq" "senr" "sent" "senv" "seo" "serf" "sers" "set" "setc" "setu" "sev" "sew" "sey" "sfg" "sfw" "shew" "shf" "shl" "shn" "shw" "shx" "sja" "sks" "slb" "slj" "slq" "slt" "smac" "smar" "smav" "smin" "smm" "smp" "smur" "spe" "spis" "sply" "spq" "spsw" "sra" "sre" "srs" "sry" "ssz" "stm" "sty" "sulf" "svo" "swd" "swp" "tad" "tan" "tbd" "tcc" "tci" "tcl" "tdl" "tet" "tgo" "thap" "thj" "thk" "tht" "thu" "tmc" "tmz" "tot" "tpn" "tpq" "tps" "tpty" "tpv" "tsp" "ttt" "vaa" "vaf" "vag" "vam" "van" "var" "vau" "vbo" "vca" "vcy" "vdb" "vei" "vex" "vff" "vfi" "vfl" "vfu" "vga" "vha" "vhr" "vit" "vmi" "vna" "vnl" "vow" "vpa" "vpb" "vpd" "vpf" "vqi" "vra" "vro" "vsa" "vsc" "vsp" "vta" "vtu" "vvl" "vvy" "wbm" "wcl" "wed" "wen" "weo" "wol" "woo" "wri" "xac" "xao" "xbo" "xbv" "xca" "xcb" "xcc" "xcf" "xci" "xcj" "xcm" "xcn" "xcp" "xct" "xcu" "xcw" "xdo" "xfh" "xfl" "xfs" "xfu" "xho" "xne" "xnm" "xom" "xoo" "xop" "xor" "xoy" "xoz" "xpo" "xtw" "yak" "yal" "yee" "yef" "yen" "yep" "yet" "yew" "yfr" "yhi" "yin" "ykr" "yma" "ypb" "ypc" "ypf" "ypi" "ypq" "ypr" "ypu" "ypy" "yrb" "ysi" "zpa" "zpl" "bthy" "smeg" "rue" "rpd" "scon" "spm" "sphu" "lio" "abq" "sut" "bwe" "smel" "scos" "spg" "rua" "smer" "swi" "soi" "rpe" "sps" "lwi" "suq" "abf" "bww" "sik" "rpt" "sph" "bths" "suz" "ati" "bmyo" "sphm" "siq" "rpx" "esi" "sud" "smd" "rmb" "sio" "sux" "rhi" "rsp" "siz" "bty" "spj" "stax" "ahu" "nwi" "eat" "bmyc" "spk" "sphi" "azt" "ean" "suw" "slu" "sfh" "rsh" 
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