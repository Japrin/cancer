#!/bin/bash

#export MODULESHOME=/usr/share/Modules
#. /usr/share/Modules/init/bash
export MODULEPATH="/tools/modulefiles:/workspace/zhengliangtao/05.setting/modulefiles":$MODULEPATH

module load fastp/0.23.4
module load seqtk/1.4
module load bwa/0.7.17-r1198-dirty
module load samtools/1.17
module load java/18.0.1.1
module load picard/3.0.0
module load gatk/4.2.6.1
#module load samtools/1.14
module load htslib/github
module load varScan/2.4.4
module load strelka2/2.9.10
module load strelka/1.0.14
module load bedtools/2.30.0
module load annovar/20220607



#module load CGpipeline/1.0
###module load glibc/2.12.2
#export PKG_CONFIG_PATH="/lustre1/zeminz_pkuhpc/05.setting/PKGCONFIG"

#### don't echo anything
####echo "pipeline init file was included at host:" `hostname`
export PIPELINE="/workspace/zhengliangtao/02.pipeline"

export PATH=$PIPELINE/cancer/all/:$PATH
#PATH=$PIPELINE/bin:$PATH
#PATH=$PIPELINE/bin/sjm:$PATH
#PATH=$PIPELINE/bin/gmap-gsnap:$PATH
#PATH=$PIPELINE/bin/blat:$PATH
#PATH=$PIPELINE/bin/bedtools:$PATH
#PATH=$PIPELINE/bin/vcftools:$PATH
#PATH=$PIPELINE/bin/annovar:$PATH
#PATH=$PIPELINE/bin/bamUtil:$PATH
#PATH=$PIPELINE/bin/samtools:$PATH
#PATH=$PIPELINE/bin/cnvnator:$PATH
#PATH=$PIPELINE/bin/breakdancer:$PATH
#PATH=$PIPELINE/bin/tabix:$PATH
#PATH=$PIPELINE/bin/bwa:$PATH
#PATH=$PIPELINE/bin/bowtie1:$PATH
#PATH=$PIPELINE/bin/bowtie2:$PATH
#PATH=$PIPELINE/bin/bismark:$PATH
#PATH=$PIPELINE/bin/tophat:$PATH
#PATH=$PIPELINE/bin/cufflinks:$PATH
#PATH=$PIPELINE/bin/seqtk:$PATH
#PATH=$PIPELINE/sge:$PATH
#PATH=/Share/BP/zhenglt/01.bin/RSEM:$PATH
#
#PATH=/Share/BP/zhenglt/01.bin/R/R-3.1.1/mybuild/bin:$PATH
#PATH=/Share/BP/zhenglt/01.bin/python/Python-2.7.8/mybuild/bin:$PATH
#PATH=/usr/java/latest/bin:$PATH
#
#PATH=$ROOTSYS/bin:$PATH
#PATH=/Share/BP/zhenglt/03.toolkit/mytoolkit:$PATH
#PATH=/Share/BP/zhenglt/01.bin/ucsc/linux.x86_64/:$PATH
#export PATH=$PATH
#
#export PYTHONPATH=/Share/BP/zhenglt/04.lib/python/lib/python2.7/site-packages
#
#PERL5LIB=/Share/BP/zhenglt/04.lib/perl/share/perl5:/Share/BP/zhenglt/04.lib/perl/lib64/perl5:/Share/BP/zhenglt/04.lib/perl/lib/perl5:/Share/BP/zhenglt/04.lib/perl/lib/perl5/x86_64-linux-thread-multi:${PERL5LIB:-""}
#export PERL5LIB=$PERL5LIB
#
#export POSTGRES_HOME=/Share/BP/zhenglt/01.bin/postgreSQL/postgresql-9.3.5/mybuild
#
#export R_LIBS_USER="/Share/BP/zhenglt/04.lib/R"
#export R_PROFILE="$PIPELINE/cancer/parameter/.Rprofile"
#export R_USER="$PIPELINE/cancer/parameter/.Renviron"
#
#LD_LIBRARY_PATH=/Share/BP/zhenglt/04.lib/C/lib:${LD_LIBRARY_PATH:-""}
#LD_LIBRARY_PATH=$ROOTSYS/lib/root/:${LD_LIBRARY_PATH:-""}
#LD_LIBRARY_PATH=/Share/BP/zhenglt/01.bin/meerkat/Meerkat/src/mybamtools/lib:${LD_LIBRARY_PATH:-""}
#LD_LIBRARY_PATH=/Share/BP/zhenglt/04.lib/C/boost/boost_1_47_0/mybuild/lib:${LD_LIBRARY_PATH:-""}
#LD_LIBRARY_PATH=/Share/BP/zhenglt/01.bin/bamtools/mybuild/lib/bamtools:${LD_LIBRARY_PATH:-""}
#LD_LIBRARY_PATH=/Share/BP/zhenglt/01.bin/hdf5-1.8.14/mybuild/lib:${LD_LIBRARY_PATH:-""}
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH

#export picardDIR=/WPSnew/zhenglt/01.bin/var/picard/current
#export PICARD=$picardDIR
#export GATK=/WPSnew/zhenglt/01.bin/var/gatk/current
#export gatkJAR=$GATK/GenomeAnalysisTK.jar
#export gatkKey=$GATK/20952006_zju.edu.cn.key
#export GATKKey=$gatkKey
#export varScanDIR=/WPSnew/zhenglt/01.bin/var/varScan
##export BREAKDANCER=$PIPELINE/bin/breakdancer
#export ANNOVAR=/WPSnew/zhenglt/01.bin/var/annovar/current

bundleDir=/workspace/zhengliangtao/00.database/broad/bundle/2.8/hg38
refData=$bundleDir/Homo_sapiens_assembly38.fasta
REF_ALT=$bundleDir/Homo_sapiens_assembly38.fasta.64.alt
REF=$refData
#refDataByChr=$bundleDir/byChr
#refNFile=$bundleDir/all.NBlock.bed
knownSites=$bundleDir/Homo_sapiens_assembly38.dbsnp138.vcf
HumanDB=/workspace/zhengliangtao/01.bin/ANNOVAR/annovar/humandb

#export BLASTDB="/Share/BP/zhenglt/00.database/blastDatabase"
#
#TMPDIR="/Share/BP/zhenglt/tmp"
