#!/bin/bash

iniFile="/WPS/GR/zhengliangtao/02pipeline/novo.med.seq/cancer/parameter/init_human.sh"

while getopts c: opt
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
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
        echo "Usage: $0 [-c iniFile] <output> <file1> ..."
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 2 ]
then
        echo "Usage: $0 [-c iniFile] <output> <file1> ..."
        exit 1
fi

source $iniFile

out=`cd \`dirname $1\`; pwd`/`basename $1`
shift 1

if [ $# -lt 2 ]
then 
	printf "$out:\n\t$1\n"
	ln -f $1 $out
else
	echo "$out:"
	for f in $*
	do
		f=`cd \`dirname $f\`; pwd`/`basename $f`
		printf "\t$f";
	done
	printf "\n";
	printf "not ready......\n";
fi
