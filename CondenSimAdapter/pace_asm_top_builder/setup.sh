#!/bin/bash
#usage: path to gromacs topology directory (e.g., "/usr/local/gromacs/share/gromacs/top/")


# copy FF file
#mkdir "$1"pace-asm.ff
#cp pace-asm.ff/* "$1"pace-asm.ff/

# compile genPACEpair
cd scripts/genpair/
make
cd ../../


# generate run script
SPATH="$PWD/scripts/"
#print $SPATH
sed -e "s|TO_BE_REPLACED|$SPATH|g" -e "s|MDP_DIR|$PWD|g" "$SPATH"prepare_temp.sh > prepare.sh || exit


