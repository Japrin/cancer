#!/bin/bash -eu

echo "*** Performing Indel Analysis using the Dindel model with the GATK Unified Genotyper ***"

if [ $# -lt 2 ]
then 
	echo "Usage: $0 <output> <bam>..."
	exit 1
fi

source analysis_20130516/Z001/somatic/Tumor/varScan/Z001.Tumor.varScan.snp.vcf

o=`cd \`dirname $1\`; pwd`/`basename $1`
shift
o_dir=`dirname $o`

f=''
for i in $*
do
        f="$f -I `cd \`dirname $i\`; pwd`/`basename $i`"
done

optL=''
if [ $# -eq 1 ]
then
        if [[ "$1" =~ .*chr[^.]*\..* ]]
        then
                chr=`echo "$1" | sed 's/.*\(chr[^\.]*\)\..*/\1/'`
		optL="-L $chr"
        fi
fi

echo ">>> Running the unified genotyper for indel calling"
java -Xms5g -Xmx5g -Djava.io.tmpdir=$o_dir -jar $GATK/GenomeAnalysisTK.jar \
	-T UnifiedGenotyper \
	$f \
	-R $REF \
	-o $o $optL \
	-A AlleleBalance \
	-A DepthOfCoverage \
	-A MappingQualityZero \
	-baq CALCULATE_AS_NECESSARY \
 	-stand_call_conf 30.0 \
	-stand_emit_conf 10.0 \
	-glm INDEL \
	-et NO_ET \
	-K $GATKKey

echo "*** Finished Indel Analysis using the GATK Unified Genotyper ***"
