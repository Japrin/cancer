#!/bin/bash

iniFile="/PROJ/GR/share/medinfo.02pipeline/cancer/parameter/init_human.sh"
TR=""

while getopts c:r: opt
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
	r)	
		if [ -f $OPTARG ]
		then
			TR="$OPTARG"
		else
			echo "WARNING: invalid target file ($OPTARG), no target will be used"
		fi
		;;

	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-r target region] <in> <out>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 2 ]
then 
	echo "Usage: $0 [-c iniFile] [-r target region] <in> <out>"
	exit 1
fi

source $iniFile


echo begin at: `date` PID: $$

ifile=$1
ofile=$2

oDir=`dirname $ofile`
mkdir -p $oDir

source $ifile

toRM=""

#### data production
##$dataProductTable
dataProductStat.sh $dataProduction_list $oDir/dataProduction.txt
#### mapping and coverage
collectMappingCov.sh ${stat_list} $oDir/mapping_coverage.txt

#### depth distribution
collectDepthPlot.pl -i=0 -j=2 $depth_list   > $oDir/depth.all.for.plot
collectDepthPlot.pl -i=0 -j=1 -t $cumu_list > $oDir/cumu.all.for.plot
Rscript $PIPELINE/all/depth.plot.R $oDir/depth.all.for.plot $oDir/cumu.all.for.plot $oDir/depth.png
#toRM="$toRM $oDir/depth.all.for.plot $oDir/cumu.all.for.plot"

#### coverage by chromosome
collectCovByChr.sh $covByChr_list $oDir/covByChr.png
#toRM="$toRM $oDir/covByChr.png.1.txt $oDir/covByChr.png.2.txt"


#### snp stat
snpList=$oDir/snp.list
awk '/sample_GATK\.snp\./' $var_sample_GATK_list > $snpList
collectVarStat.sh $snpList $oDir/snp.stat.txt
awk 'NR==1 || !/^sampleID/' $oDir/snp.stat.txt.genome.txt $oDir/snp.stat.txt.exome.txt > $oDir/snp.stat.txt
#rm $oDir/snp.stat.txt.genome.txt $oDir/snp.stat.txt.exome.txt $snpList

#### indel stat
indelList=$oDir/indel.list
awk '/sample_GATK\.indel\./' $var_sample_GATK_list > $indelList
collectVarStat.sh $indelList $oDir/indel.stat.txt indel
awk 'NR==1 || !/^sampleID/' $oDir/indel.stat.txt.genome.txt $oDir/indel.stat.txt.exome.txt > $oDir/indel.stat.txt
rm $oDir/indel.stat.txt.genome.txt $oDir/indel.stat.txt.exome.txt $indelList

#### generate 
printf "" > $oDir/gen.in
printf "data_production=dataProduction.txt\n" >> $oDir/gen.in
printf "mapping_stat=mapping_coverage.txt\n" >> $oDir/gen.in
printf "snp_stat=snp.stat.txt\n" >> $oDir/gen.in
printf "indel_stat=indel.stat.txt\n" >> $oDir/gen.in

aaDir=`pwd`
cd $oDir 
cp $PIPELINE/report/001.jpg .
cp $PIPELINE/report/report_template.bib .
cp $PIPELINE/report/background.pdf .
report_gen.pl -l gen.in -t $PIPELINE/report/report_template.tex report.tex
#/WPS/GR/zhengliangtao/01bin/latext/2012/bin/x86_64-linux/xelatex -interaction=nonstopmode report.tex
#/WPS/GR/zhengliangtao/01bin/latext/2012/bin/x86_64-linux/bibtex report
#/WPS/GR/zhengliangtao/01bin/latext/2012/bin/x86_64-linux/xelatex -interaction=nonstopmode report.tex
#/WPS/GR/zhengliangtao/01bin/latext/2012/bin/x86_64-linux/xelatex -interaction=nonstopmode report.tex
#rm 001.jpg report_template.bib background.pdf
#cd $aaDir
#mv $oDir/report.pdf $ofile


rm -f $toRM


echo end at: `date`
