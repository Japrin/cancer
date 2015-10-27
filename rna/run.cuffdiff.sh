#!/bin/bash


iniFile="/WPS/BP/zhenglt/02.pipeline/health.02pipeline/cancer/parameter/init_human.sh"
while getopts c:t:f:ab opt
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
	        echo "Usage: $0 [-c iniFile]  <outDir> <listFile>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 2 ]
then 
	echo "Usage: $0 [-c iniFile]  <outDir> <listFile>"
	exit 1
fi



source $iniFile

genome="/WPS/BP/zhenglt/00.database/broad/bundle/2.8/b37/bowtie2/human_g1k_v37_decoy"
annGTF="/WPS/BP/zhenglt/00.database/ensemble/release75/homo_sapiens/gtf/Homo_sapiens.GRCh37.75.gtf"
transcriptome="/WPS/BP/zhenglt/00.database/broad/bundle/2.8/b37/bowtie2/human_g1k_v37_decoy.transcriptome/Homo_sapiens.GRCh37.75"
IGVGenome="/WPS/BP/zhenglt/00.database/broad/bundle/2.8/b37/human_g1k_b37_decoy5.igv.genome"

outDir=$1
listFile=$2
#lable1=$3
#lable2=$4

echo begin at: `date`

mkdir -p $outDir

####cxbList1=`awk '/'$lable1'/{print $3}' $listFile | xargs | tr ' ' ','`
####cxbList2=`awk '/'$lable2'/{print $3}' $listFile | xargs | tr ' ' ','`


#echo "cuffdiff --no-update-check -p 8 -o $outDir --use-sample-sheet $annGTF $listFile"
#cuffdiff --no-update-check -p 8 -o $outDir --use-sample-sheet $annGTF $listFile

echo "cuffnorm --no-update-check -p 8 -o $outDir --use-sample-sheet $annGTF $listFile"
cuffnorm --no-update-check -p 8 -o $outDir --use-sample-sheet $annGTF $listFile

#awk 'NR==1||$NF=="yes"' $outDir/gene_exp.diff > $outDir/gene_exp.sig.diff
#awk 'NR==1||$NF=="yes"' $outDir/isoform_exp.diff > $outDir/isoform_exp.sig.diff
      
echo end at: `date`
