/*
 * =====================================================================================
 *
 *       Filename:  chr22.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/15/2012 11:11:42 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
{
	char * rootFile="/WPS/GR/zhengliangtao/work/CBTd/analysis/SV/CNVnator/bin100/CBTd.root";
	int BIN_SIZE=100;
	int nchr=24;
	char * CHR_LIST[nchr]={"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};

	TFile tf=TFile(rootFile);
	TDirectory * pBinDir=tf.Get("bin_100;1");
	
	for(int j=0;j<nchr;j++)
	{
		char * CHR=CHR_LIST[j];
		char hisName[256];
		sprintf(hisName,"his_rd_p_%s_100_partition_GC_merge;1",CHR);

		TH1D * pHis=pBinDir->Get(hisName);
		int nbin=pHis->GetNbinsX();
		for(int i=1;i<nbin;i++)
		{
				printf("%s\t%d\t%d\t%d\t%d\n",CHR,i,i*BIN_SIZE+1,(i+1)*BIN_SIZE,pHis->GetBinContent(i));
		}
	}
}
