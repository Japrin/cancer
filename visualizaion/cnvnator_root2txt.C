/*
 * =====================================================================================
 *
 *       Filename:  macro.cnvnator.root2txt.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/15/2012 11:11:42 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  zhengliangtao (), 
 *   Organization:  
 *
 * =====================================================================================
 */


int cnvnator_root2txt(char * rootFile,int BIN_SIZE)
{
	//char * rootFile="/PROJ/GR/HUMAN/shidao.cancer.wgs/zhengliangtao/analysis/cnvnator/out/CS-1/CS-1.cnvnator.root";
	//int BIN_SIZE=100;
	const int nchr=24;
	char * CHR_LIST[nchr]={"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};

	char * sBinDir[256];
	sprintf(sBinDir,"bin_%d;1",BIN_SIZE);
	const char* csBinDir=sBinDir;
	TFile tf=TFile(rootFile);
	TDirectory * pBinDir=tf.Get(csBinDir);
	//TDirectory * pBinDir=tf.Get("bin_100;1");
	
	for(int j=0;j<nchr;j++)
	{
		char * CHR=CHR_LIST[j];
		char hisMerge[256];
		sprintf(hisMerge,"his_rd_p_%s_%d_partition_GC_merge;1",CHR,BIN_SIZE);
		char hisPartition[256];
		sprintf(hisPartition,"his_rd_p_%s_%d_partition_GC;1",CHR,BIN_SIZE);
		char hisGC[256];
		sprintf(hisGC,"his_rd_p_%s_%d_GC_l1;1",CHR,BIN_SIZE);
		char hisGC2[256];
		sprintf(hisGC2,"his_rd_p_%s_%d_GC_l2;1",CHR,BIN_SIZE);
		char hisGC3[256];
		sprintf(hisGC3,"his_rd_p_%s_%d_GC_l3;1",CHR,BIN_SIZE);
		char hisRaw[256];
		sprintf(hisRaw,"his_rd_p_%s_%d;1",CHR,BIN_SIZE);

		TH1D * pHisMerge=pBinDir->Get(hisMerge);
		int nbinMerge=pHisMerge->GetNbinsX();
		TH1D * pHisPartition=pBinDir->Get(hisPartition);
		int nbinPartition=pHisPartition->GetNbinsX();
		TH1D * pHisGC=pBinDir->Get(hisGC);
		int nbinGC=pHisGC->GetNbinsX();
		TH1D * pHisGC2=pBinDir->Get(hisGC2);
		int nbinGC2=pHisGC2->GetNbinsX();
		TH1D * pHisGC3=pBinDir->Get(hisGC3);
		int nbinGC3=pHisGC3->GetNbinsX();
		TH1D * pHisRaw=pBinDir->Get(hisRaw);
		int nbinRaw=pHisRaw->GetNbinsX();
		printf("#nbinMerge: %d\tnbinPartition: %d\tnbinGC: %d\tnbinGC2: %d\tnbinGC3: %d\tnbinRaw: %d\n",nbinMerge,nbinPartition,nbinGC,nbinGC2,nbinGC3,nbinRaw);
		printf("#chrom\tbeg\tend\ti\tRD_merge\tRD_partition\tRD_GC1\tRD_GC2\tRD_GC3\tRD_raw\n");
		for(int i=1;i<nbinMerge;i++)
		{
				//printf("%s\t%d\t%d\t%d\t%4.2f\t%4.2f\n",CHR,i*BIN_SIZE+1,(i+1)*BIN_SIZE,i,pHisMerge->GetBinContent(i),pHisGC->GetBinContent(i));
				printf("%s\t%d\t%d\t%d\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\n",CHR,i*BIN_SIZE+1,(i+1)*BIN_SIZE,i,pHisMerge->GetBinContent(i),pHisPartition->GetBinContent(i),pHisGC->GetBinContent(i),pHisGC2->GetBinContent(i),pHisGC3->GetBinContent(i),pHisRaw->GetBinContent(i));
		}
	}
}
