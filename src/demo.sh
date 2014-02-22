#!/bin/bash

set -e

export PATH=$PATH:/usr/bin:/bin

# check number of input arguments
if [ $# != "5" ]; then
	echo -e "usage:\n\t `basename $0` inpattern FIRST LAST NSEEDS outprefix"
	#                                 1         2     3    4      5
	exit 1
fi

# assign input arguments
IN_PATTERN=$1
FIRST=$2
LAST=$3
NSEEDS=$4
OUT_PREFIX=$5

echo "IN_PATTERN=$1"
echo "FIRST=$2"
echo "LAST=$3"
echo "NSEEDS=$4"
echo "OUT_PREFIX=$5"

INFILES=""
for i in `seq $FIRST $LAST`; do
	F=`printf "$IN_PATTERN" $i`
	INFILES="$INFILES $F"
done
imcombine avg $INFILES -o ${OUT_PREFIX}iavg.png &
imcombine weisz $INFILES -o ${OUT_PREFIX}imed.png &

if [ "$NSEEDS" -lt 1 ]; then NSEEDS=1 ; fi
if [ "$NSEEDS" -gt 9 ]; then NSEEDS=9 ; fi
# set-up num-threads for a 32 core machine
NT=0
case "$NSEEDS" in
"1") NT=30 ;;
"2") NT=15 ;;
"3") NT=10 ;;
"4") NT=7 ;;
"5") NT=6 ;;
"6") NT=5 ;;
"7") NT=4 ;;
"8") NT=4 ;;
"9") NT=3 ;;
esac

NPROC=`nproc`
if [ "$NT" -gt "$NPROC" ]; then NT=$NPROC ; fi
if [ "$FORCE_NT" -gt "0" ]; then NT=$FORCE_NT ; fi

echo "NT=$NT"


CENTROIDS=""
for i in `seq 0 $[NSEEDS-1]`; do
	SEED=`printf "$IN_PATTERN" $[FIRST+i]`
	CENTROID=${OUT_PREFIX}c$i.png
	OMP_NUM_THREADS=$NT centroid $IN_PATTERN $FIRST $LAST $SEED $CENTROID &
	CENTROIDS="$CENTROIDS $CENTROID"
done
wait
echo "CENTROIDS=\"$CENTROIDS\""
imcombine avg $CENTROIDS -o ${OUT_PREFIX}cavg.png
imcombine weisz $CENTROIDS -o ${OUT_PREFIX}cmed.png
