#!/bin/bash

while getopts t: opt
do
	case $opt in 
	t)
		optT=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
        echo "Usage: $0 <outDir> <fq1> <fq2> <sampleID>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
    echo "Usage: $0 <outDir> <fq1> <fq2> <sampleID>"
	exit 1
fi

outDir=$1
fq1=$2
fq2=$3
sampleID=$4
mkdir -p $outDir

export HLAminerDir=/Share/BP/zhenglt/01.bin/hla/hlaminer/HLAminer_v1.3.1
ncbiBlastConfigFile="/Share/BP/zhenglt/02.pipeline/cancer/MHC/HLAminer/ncbiBlastConfig2-2-22.txt"

fof_file=$outDir/patient.fof
printf "$fq1\n$fq2\n" > $fof_file

echo begin at: `date`
cd $outDir
###Run TASR
echo "Running TASR..."
#TASR Default is -k 15 for recruiting reads. You may increase k, as long as k < L/2 where L is the minimum shotgun read length
$HLAminerDir/bin/TASR -f $fof_file -m 20 -k 20 -s $HLAminerDir/database/HLA-I_II_GEN.fasta -i 1 -b TASRhla -w 1
###Restrict 200nt+ contigs
cat TASRhla.contigs |perl -ne 'if(/size(\d+)/){if($1>=200){$flag=1;print;}else{$flag=0;}}else{print if($flag);}' > TASRhla200.contigs
###Create a [NCBI] blastable database
echo "Formatting blastable database..."
$HLAminerDir/bin/formatdb -p F -i TASRhla200.contigs
###Align contigs against database
echo "Aligning TASR contigs to HLA references..."
$HLAminerDir/bin/parseXMLblast.pl -c $ncbiBlastConfigFile -d $HLAminerDir/database/HLA-I_II_GEN.fasta -i TASRhla200.contigs -o 0 -a 1 > tig_vs_hla-ncbi.coord
###Align HLA references to contigs
echo "Aligning HLA references to TASR contigs (go have a coffee, it may take a while)..."
$HLAminerDir/bin/parseXMLblast.pl -c $ncbiBlastConfigFile -i $HLAminerDir/database/HLA-I_II_GEN.fasta -d TASRhla200.contigs -o 0 > hla_vs_tig-ncbi.coord
###Predict HLA alleles
echo "Predicting HLA alleles..."
$HLAminerDir/bin/HLAminer.pl -b tig_vs_hla-ncbi.coord -r hla_vs_tig-ncbi.coord -c TASRhla200.contigs -h $HLAminerDir/database/HLA-I_II_GEN.fasta -p $HLAminerDir/database/hla_nom_p.txt

echo end at: `date`
