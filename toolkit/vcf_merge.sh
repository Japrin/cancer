#!/bin/bash -eu

iniFile="`dirname $0`/../parameter/init_human.sh"

while getopts c: opt
do
	case $opt in 
	c)	
		if [ -f $OPTARG ]
		then
			iniFile="$OPTARG"
		else
			echo "WARNING: invalid configure file ($OPTARG), default will be used"
		fi
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] <output> <VCF file to merge>..."
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 2 ]
then 
	echo "Usage: $0 [-c iniFile] <output> <VCF file to merge>..."
	exit 1
fi
source $iniFile 
output=`cd \`dirname $1\`; pwd`/`basename $1`

shift

inputs=''
for i in $*
do
	i=`cd \`dirname $i\`; pwd`/`basename $i`
	echo ">> Indexing VCF: $i"
	if [ ! -f "$i.gz" ]
	then
		bgzip -f $i
	fi
	tabix -f -p vcf $i.gz
	inputs="$inputs $i.gz"
done

echo begin at: `date`
echo "*** Performing VCF Merging ***"
echo ">> Merging VCFs into $output"

if [ ! -f $output ]
then
	vcf-concat $inputs > $output
	#vcf-merge $inputs > $output
fi
