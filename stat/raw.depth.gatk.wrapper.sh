#!/bin/bash

echo "***  ***"

iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"

while getopts c:w:g: opt
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
		echo "Usage: $0 [-c iniFile] <prefix>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 1 ]
then 
	echo "Usage: $0 [-c iniFile] <prefix>"
	exit 1
fi

echo begin at: `date`

prefix=$1

source $iniFile

printf "*** hostname: %s ***\n" `hostname`

sed -e  '2,$s/[:-]/\t/g' -e '1,1s/Target/chr\tbeg\tend/' $prefix.sample_interval_summary > $prefix.sample_interval_summary.forPlot
Rscript $PIPELINE/cancer/stat/raw.depth.gatk.R $prefix.sample_interval_summary.forPlot

hfile=$prefix.sample_interval_statistics
printf "options(bitmapType='cairo')\na<-t(read.table(\"$hfile\",header=T))\npng(\"$hfile.png\",width=800,height=600)\nplot(x=0:(dim(a)[1]-1),y=a[,1],type='p',pch=20,col='darkred',xlab='Cumulative Depth',ylab='Count',main='Distribution of per region')\ndev.off()\n" | R --vanilla --slave

hfile=$prefix.sample_statistics
printf "options(bitmapType='cairo')\na<-t(read.table(\"$hfile\",header=T))\npng(\"$hfile.png\",width=800,height=600)\nplot(x=0:(dim(a)[1]-1),y=a[,1],type='p',pch=20,col='darkred',xlab='Depth',ylab='Count',main='Distribution of per base depth')\ndev.off()\n" | R --vanilla --slave


echo end at: `date`
