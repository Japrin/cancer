#!/bin/bash

iniFile="/Share/BP/zhenglt/02.pipeline/cancer/parameter/init_human.sh"
optV=""
while getopts c:v: opt
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
    v)
        if [[ $OPTARG != "" ]];then
            optV="-v $OPTARG"
        fi
        ;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-c iniFile] [-v virus file] <sampleID> [<fq1> <fq2>|bam] <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-v virus file] <sampleID> [<fq1> <fq2>|bam] <outDir>"
	exit 1
fi

echo begin at: `date`

source $iniFile
module unload samtools/1.3.1
module load samtools/0.1.19
#module list

configTemplateFile=/Share/BP/zhenglt/01.bin/virusIntegration/VirusFinder/VirusFinder2_with_VERSE/VirusFinder2.0/template-config.RNA.b37.txt
VirusFinderBin=/Share/BP/zhenglt/01.bin/virusIntegration/VirusFinder/VirusFinder2_with_VERSE/VirusFinder2.0/VirusFinder.pl

if [ $# -eq 3 ];then
    sampleID=$1
    bamfile=$2
    outDir=$3
    mkdir -p $outDir
    configFile=$outDir/$sampleID.config.txt
    perl -ane 'chomp;
           if(/^#alignment_file/){print "alignment_file = '$bamfile'\n";}
           else{ print "$_\n";} ' $configTemplateFile > $configFile
elif [ $# -eq 4 ];then
    sampleID=$1
    fq1=$2
    fq2=$3
    outDir=$4
    mkdir -p $outDir
    configFile=$outDir/$sampleID.config.txt
    perl -ane 'chomp;
           if(/^#fastq1/){print "fastq1         = '$fq1'\n";}
           elsif(/^#fastq2/){print "fastq2         = '$fq2'\n";}
           else{ print "$_\n";} ' $configTemplateFile > $configFile
fi
echo perl $VirusFinderBin $optV -c $configFile -o $outDir
perl $VirusFinderBin $optV -c $configFile -o $outDir >$outDir/$sampleID.log 2>&1

echo end at: `date`
