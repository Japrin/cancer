#!/bin/bash -eu

echo "*** Performing Variant Filtration ***"

if [ $# -lt 1 ]
then 
	echo "Usage: $0 <vcf>"
	exit 1
fi

source /WPS/GR/zhengliangtao/02pipeline/novo.med.seq/cancer/parameter/init_human.sh

f=`cd \`dirname $1\`; pwd`/`basename $1`
o=${f/.vcf/}.filtered.vcf
filter="(AB ?: 0) > 0.75 || QUAL < 50.0 || DP > 500 || SB > -0.1 || MQ0>=4"
o_dir=`dirname $o`

echo ">> Filtering variants using filter: $filter"
java -Xms5g -Xmx5g -Djava.io.tmpdir=$o_dir -jar $GATK/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R $REF \
	-o $o \
	-V $f \
	--clusterWindowSize 10 \
	--filterExpression "$filter" \
	--filterName "StandardFilter" \
	-et NO_ET \
	-K $GATKKey

mv $o $f
rm -f $o.*

echo "*** Finished Variant Filtration ***"
