 
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
  "cyo" "dak" "dal" "dar" "das" "dat" "dba" "ddi" "def" "del" "deo" "dfa" "dfl" "dgg" "dgt" "dha" "dhk" "dhy" "dia" "dma" "doa" "dpg" "dpi" "dpp" "dps" "drt" "dsa" "dsf" "dto" "dtr" "dts" "dzi" "eaa" "eae" "eal" "ear" "eba" "ech" "echg" "echj" "echl" "echs" "echv" "echw" "ecoi" "ecoj" "ecor" "ecos" "edi" "eex" "efe" "egl" "egr" "ehe" "ehi" "eic" "ein" "eiv" "ema" "eme" "emr" "epa" "erc" "erg" "ero" "eru" "erw" "fam" "fcy" "fso" "fta" "ftf" "ftg" "fth" "fti" "ftl" "fto" "ftq" "ftr" "fts" "ftt" "ftv" "gan" "gho" "gmx" "gra" "gsj" "gsl" "haa" "haf" "hbe" "hhh" "hhs" "hmg" "hol" "hro" "hvn" "idt" "ili" "ilo" "ipi" "jeu" "kaf" "kak" "kbl" "kbt" "kci" "kct" "kde" "kga" "kgy" "kla" "kll" "klw" "kma" "kmi" "kmx" "koc" "koe" "kok" "kox" "kpa" "kpb" "kpc" "kpe" "kpg" "kph" "kpi" "kpj" "kpk" "kpm" "kpne" "kpnk" "kpnu" "kpo" "kpq" "kpr" "kpt" "kpw" "kpy" "kpz" "kqu" "kqv" "kuy" "kva" "kvq" "lag" "lcd" "lel" "les" "lfa" "lgi" "lgt" "lha" "lhk" "lib" "lip" "llg" "llo" "loa" "lpc" "lpe" "lpf" "lpm" "lpp" "lsd" "lsh" "lth" "mah" "mbac" "mbr" "mbur" "mca" "mdh" "mdn" "meh" "mei" "mep" "metl" "metr" "meu" "mfa" "mgl" "mko" "mlo" "mmai" "mmas" "mmk" "mmn" "mpsy" "mrt" "msd" "mtm" "mxa" "myi" "naf" "nai" "nba" "nce" "ncr" "ncs" "ndi" "nei" "nek" "nel" "ngd" "ngk" "ngo" "nhm" "nma" "nmc" "nmj" "nmn" "nmq" "nmt" "nmw" "nmx" "nri" "nse" "nsi" "nte" "nve" "nwe" "odi" "ott" "our" "ovi" "paca" "pacr" "paei" "paf" "pag" "pagc" "pagg" "palh" "palk" "paln" "paly" "pan" "panp" "panr" "pans" "pant" "pap" "part" "pat" "pau" "paur" "pay" "pbb" "pbc" "pbw" "pcan" "pcd" "pch" "pck" "pcp" "pcr" "pcy" "pcz" "pdam" "pdr" "pds" "pea" "pfa" "pfd" "pfh" "pfk" "pfz" "pgb" "pgr" "pgz" "pha" "phei" "pia" "pic" "pif" "pig" "pkc" "pkn" "pkr" "pkz" "plb" "plc" "pld" "ple" "pli" "plo" "plq" "plu" "plum" "ply" "pmai" "pmib" "pnc" "pol" "pos" "ppa" "pphe" "pprf" "ppsy" "ppuu" "ppv" "prj" "prq" "prw" "pse" "psec" "psed" "psem" "psen" "pseo" "pset" "psi" "psil" "psk" "pso" "psos" "pspg" "pspo" "pspw" "psta" "pstw" "psv" "psw" "psx" "psyg" "ptc" "pti" 
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