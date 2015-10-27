#!/bin/bash

opt_tad=0.1
opt_g=1
opt_v=1.5

while getopts c:t:g:v: opt
do
	case $opt in 
	c)	
		if [ -f $OPTARG ]
		then
			iniFile="$OPTARG"
		else
			echo "WARNING: invalid reference file ($OPTARG), default will be used"
		fi
		;;
	t)
		opt_tad=$OPTARG
		;;
	g)
		opt_g=$OPTARG
		;;
	v)
		opt_v=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-t amplification and deletion threshold, defautl 0.1] [-g geneGistic, default 1] [-v cap, default 1.5] <seg.file> <marker.file> <out.dir> <analysis.id>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-t amplification and deletion threshold, defautl 0.1] [-g geneGistic, default 1] [-v cap, default 1.5] <seg.file> <marker.file> <out.dir> <analysis.id>"
	exit 1
fi

# --- SET UP ENVIRONMENT VARIABLES ---
# presumed location of MATLAB Component Runtime (MCR) v7.14
# if the MCR is in a different location, edit the line below
GISTICDir="/Share/BP/zhenglt/01.bin/gistic/GISTIC_2_0_22"
mcr_root="/Share/BP/zhenglt/01.bin/gistic/GISTIC_2_0_22/MATLAB_Component_Runtime"

LD_LIBRARY_PATH=$mcr_root/v714/runtime/glnxa64:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=$mcr_root/v714/sys/os/glnxa64:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=$mcr_root/v714/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=$mcr_root/v714/sys/java/jre/glnxa64/jre/lib/amd64/server:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=$mcr_root/v714/sys/java/jre/glnxa64/jre/lib/amd64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export XAPPLRESDIR=$mcr_root/v714/X11/app-defaults
# (these may be set permanently by copying the above lines into your login script)

# --- RUN GISTIC 2.0 ---
segfile=$1
markersfile=$2
outDir=$3
analysis_id=$4
refgenefile=$GISTICDir/refgenefiles/hg19.mat

echo begin at: `date`

mkdir -p $outDir

$GISTICDir/gp_gistic2_from_seg -b $outDir -seg $segfile -mk $markersfile -refgene $refgenefile \
	-genegistic $opt_g -smallmem 1 -broad 1 -brlen 0.5 -conf 0.99 -armpeel 1 -savegene 1 -gcm mean -fname $analysis_id \
	-ta $opt_tad -td $opt_tad -cap $opt_v

echo end at: `date`
