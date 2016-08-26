#!/usr/bin/env python
import h5py
import sys
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import cm
from scLVM import scLVM
import limix.modules.panama as PANAMA
import limix.modules.varianceDecomposition as VAR
import limix.modules.qtl as QTL
from scLVM.utils.misc import *
from scLVM.utils.barplot import *
import scipy as SP
import limix
import limix.modules.varianceDecomposition as VAR
import pdb
import traceback
import pylab as PL
import argparse
import os
import glob
import numpy as np

if __name__ == '__main__':

    version = "1.0"
    parser = argparse.ArgumentParser(description='Run scLVM analysis on data in hdf5 format')
    parser.add_argument("-o", "--output", metavar="STR", type=str, dest="outDir",
                        default=".", help="output directory [default: %(default)s]")
    parser.add_argument("-i", "--input", metavar="STR", type=str, dest="inFile",
                        required=True, help="input hdf5 file [required]")
    parser.add_argument("-p", "--panama", metavar="STR", type=str, dest="panama",
                        required=True, help="cached file, cell-cycle covariance matrix in hdf5 file format [required]")
    parser.add_argument("-N", "--NJobs", metavar="INT", type=int, dest="NJobs",
                        default=1, help="number of jobs to split into [default: %(default)s]")
    parser.add_argument("-j", "--jJob", metavar="INT", type=int, dest="jJob",
                        default=1, help="number of jobs to split into [default: %(default)s]")
    parser.add_argument("-s", "--step", metavar="STR", type=str, dest="step",
                        default="fit", help="scLVM step, possible steps:(fit, decomp, collect) [default: %(default)s]")
    parser.add_argument("-a", "--sample", metavar="STR", type=str, dest="sample",
                        default="SAMPLE", help="sample id [default: %(default)s]")
    args = parser.parse_args()
    NJobs = args.NJobs                          # number of jobs to split to
    jJob = args.jJob - 1                        # 0-based in program
    #Where to save results
    out_name = 'correlation_test'
    out_dir = os.path.join(args.outDir,out_name)
    run_dir = os.path.join(out_dir,'runsTcells')
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)
    #load data
    f = h5py.File(args.inFile,'r')
    Y = f['Y'][:]
    tech_noise = f['tech_noise'][:]             # technical noise
    genes_het_bool=f['genes_het_bool'][:]       # index of heterogeneous(??!??) genes
    geneID = f['geneID'][:]                     # gene names
    cellcyclegenes_filter = SP.unique(f['cellcyclegenes_filter'][:].ravel() - 1)    # idx of cell cycle genes
    cellcyclegenes_filterCB600 = f['cellcyclegenes_filterCB'][:].ravel() - 1        # idxof cell cycle genes ...

    # filter cell cycle genes
    idx_cell_cycle = SP.union1d(cellcyclegenes_filter, cellcyclegenes_filterCB600)
    Ymean2 = Y.mean(0)**2 > 0
    idx_cell_cycle_noise_filtered = SP.intersect1d(idx_cell_cycle,SP.array(SP.where(Ymean2.ravel() > 0)))
    Ycc = Y[:, idx_cell_cycle_noise_filtered]
    ###pdb.set_trace()
    sclvm = scLVM(Y)
    if args.step == "fit":
        # Fit GPLVM to data
        #k = 80
        k = min(Ycc.shape[0],80)
        X_ARD,Kcc_ARD,varGPLVM_ARD = sclvm.fitGPLVM(idx=idx_cell_cycle_noise_filtered,k=k,
                                                    out_dir=os.path.join(args.outDir, 'cache'),
                                                    file_name=args.panama, recalc=True, use_ard=True)
        # Plot variance contributions from ARD
        plt = PL.subplot(1, 1, 1)
        PL.title('Variance explained by latent factors')
        PL.scatter(SP.arange(k)+1,varGPLVM_ARD['X_ARD'])
        PL.xlim([0,k+1])
        PL.xlabel('# Factor')
        PL.ylabel('Variance explained')
        PL.savefig(os.path.join(args.outDir, "variance-explained.scree.png"))

        k = 1
        X, Kcc, varGPLVM = sclvm.fitGPLVM(idx=idx_cell_cycle_noise_filtered, k=k,
                                      out_dir=os.path.join(args.outDir, 'cache'),
                                      file_name=args.panama, recalc=True, use_ard=False)
        # Plot inferred similarity matrix
        plt = PL.subplot(1, 1, 1)
        PL.title('Similarity matrix based on cell cycle')
        PL.imshow(Kcc, cmap=cm.RdBu, vmin=-3, vmax=+3, interpolation='None')
        PL.colorbar()
        plt.set_xticks([])
        plt.set_yticks([])
        PL.xlabel('cells')
        PL.ylabel('cells')
        PL.savefig(os.path.join(args.outDir, "cell-cell.similarity.png"))
    elif args.step == "decomp":
        ##if NJobs == jJob + 1:
        ##    exit("The last slice is discarded")
        ## just load covariance matrix caculated in "fit" step
        k = 1
        X, Kcc, varGPLVM = sclvm.fitGPLVM(idx=idx_cell_cycle_noise_filtered, k=k,
                                      out_dir=os.path.join(args.outDir, 'cache'),
                                      file_name=args.panama, recalc=False, use_ard=False)
        # load relevant dataset for analysis
        #genes_het = SP.array(SP.where(genes_het_bool[:].ravel() == 1))
        # considers only heterogeneous genes
        Ihet = genes_het_bool == 1
        Y = Y[:, Ihet]
        tech_noise = tech_noise[Ihet]
        geneID = geneID[Ihet]
        # split across genes
        #pdb.set_trace()
        ####Iy = SP.array(SP.linspace(0, Y.shape[1], NJobs), dtype='int')
        #### the last one will be discarded
        Iy = SP.array(SP.linspace(0, Y.shape[1], NJobs+1), dtype='int')
        i0 = Iy[jJob]
        i1 = Iy[jJob+1]
        ###i1 = Iy[min(jJob+1,NJobs-1)]
        ###i1 = max(i0, i1+1)
        #
        sclvm = scLVM(Y, geneID=geneID, tech_noise=tech_noise)
        # fit the model from i0 to i1
        #pdb.set_trace()
        #try:
        sclvm.varianceDecomposition(K=Kcc, i0=i0, i1=i1,verbose=False)
        #except:
        #    type, value, tb = sys.exc_info()
        #    traceback.print_exc()
        #    pdb.post_mortem(tb)
        #pdb.set_trace()
        # get variance components
        var, var_info = sclvm.getVarianceComponents(normalize=True)
        print(var_info['col_header'])
        print(var[0:5, ])
        # pdb.set_trace()
        # var_filtered = var[var_info['conv']] # filter out genes for which vd has not converged
        # get corrected expression levels
        Ycorr = sclvm.getCorrectedExpression()

        # fit lmm without correction
        pv0, beta0, info0 = sclvm.fitLMM(K=None, i0=i0, i1=i1, verbose=True)
        # fit lmm with correction
        pv, beta, info = sclvm.fitLMM(K=Kcc, i0=i0, i1=i1, verbose=True)
        # create outfile
        out_file = os.path.join(run_dir, 'job_%03d_%03d.hdf5' % (jJob, NJobs))
        fout = h5py.File(out_file, 'w')
        count = 0
        for i in xrange(i0, i1):
            gene_id = 'gene_%d' % (i)
            out_group = fout.create_group(gene_id)
            RV = {}
            RV['Ycorr'] = Ycorr[:, count]
            RV['pv0'] = pv0[count, :]
            RV['pv'] = pv[count, :]
            RV['beta'] = beta[count, :]
            RV['beta0'] = beta0[count, :]
            RV['vars'] = var[count, :]
            RV['varsnorm'] = var[count, :]
            # RV['vars_info'] = var_info
            RV['is_converged'] = SP.repeat(var_info['conv'][count] * 1, 7)
            ##RV['is_converged'] = SP.repeat(var_info['conv'][count] * 1, 5)
            dumpDictHdf5(RV, out_group)
            count += 1
        fout.close()
    elif args.step == "collect":
             # load data
            correlations=1    # collect data from correlation analysis
            fpa = h5py.File(os.path.join(args.outDir, 'cache', args.panama), 'r')
            cc_noise_filtered = fpa['cc_noise_filtered'][:]
            K = fpa['Kconf'][:]
            #genes_het = SP.array(SP.where(f['genes_heterogen'][:].ravel()==1)).ravel()
            genes_het = SP.array(SP.where(genes_het_bool[:].ravel() == 1)).ravel()
            # considers only heterogeneous genes
            #Ihet = genes_het_bool == 1
            Nhet = len(genes_het)
            Ncells = f['Y'][:].shape[0]
            Ngenes = f['Y'][:].shape[1]

            cc_genes=SP.ones((Ngenes,1))
            cc_genes[cc_noise_filtered]=0
            cc_genes_het=cc_genes[genes_het]

            FL = glob.glob(os.path.join(run_dir,'job*.hdf5'))

            out_file_base = os.path.join(args.outDir, args.sample + ".vars")
            out_file = out_file_base+'.hdf5'
            fout = h5py.File(out_file,'w')
            RV = {}
            if correlations ==1:
                fields = ['pv', 'pv0', 'beta', 'beta0', 'varsnorm', 'Ycorr', 'is_converged']
            else:
                fields = ['varsnorm','Ycorr','is_converged']

            for field in fields:
                RV[field] = SP.zeros([Nhet,Nhet])
            # RV['varsnorm'] = SP.zeros([Nhet,4])
            # 3 component: technical noise, cell-cycel factor, other biological
            RV['varsnorm'] = SP.zeros([Nhet,3])
            RV['Ycorr'] = SP.zeros([Nhet,Ncells])
            RV['is_converged'] = SP.zeros([Nhet,1])

            # loop through files
            for fn in FL:
                # read file
                # print "file: %s" % (fn)
                try:
                    ff = h5py.File(fn, 'r')
                    for key in ff.keys():
                        id0 = int(key.split('_')[1])
                        for field in fields:
                            if (field == 'Ycorr') or (field == 'varsnorm'):
                                RV[field][id0, :] = ff[key][field][:]
                            else:
                                if (field == 'is_converged'):
                                    RV[field][id0, 0] = ff[key][field][0, ] * 1.0
                                else:
                                    RV[field][id0, :] = ff[key][field][:]
                        pass
                    ff.close()
                except Exception:
                    #pdb.set_trace()
                    print "bad file: %s" % (fn)
                    continue
            #store
            #pdb.set_trace()
            RV['K'] = K
            RV['cc_genes_het'] = cc_genes_het
            dumpDictHdf5(RV,fout)
            fout.close()

            Ycorr = RV['Ycorr'][:]
            SP.savetxt(out_file_base + '.Ycorr.txt', Ycorr)
            SP.savetxt(out_file_base + '.isConverged.txt', RV['is_converged'])
            #plot variance decomposition
            from scLVM_cfg import *
            indsconv = RV['is_converged'].ravel() == 1
            is_nocc_conv = SP.bitwise_and(indsconv, cc_genes_het.ravel() == 1)
            vars_normConv = RV['varsnorm'][indsconv == 1, :]
            # ['hidden_0', 'biol_noise', 'tech_noise']
            # ==>
            # ['hidden_0', 'tech_noise', 'biol_noise']
            #pdb.set_trace()
            vars_normConv = vars_normConv[:, [0, 1, 2]]
            H2 = 1 - vars_normConv[:, 2]
            out_file_pdf = out_file_base+'.pdf'
            #pdb.set_trace()
            var_plot(vars_normConv, H2, CFG['var_comp_fields'][SP.array([0,1,2])],V_range=SP.linspace(0,1,11),plot_element_count=False,
                     filename=out_file_pdf, normalize=True, figsize=[6,5])
            # var_plot(vars_normConv, H2, CFG['var_comp_fields'][SP.array([0,1,2])],V_range=SP.linspace(0,1,8),plot_element_count=True,
            #         filename=out_file_pdf, normalize=True, figsize=[6,5])
            # calculate average variance components across all genes and visualize
            var_mean = vars_normConv.mean(0)
            colors = ['Green', 'MediumBlue', 'Gray']
            PL.figure(figsize=[6,5])
            pp = PL.pie(var_mean, labels=['cell cycle', 'biology', 'tech noise'],
                      autopct='%1.1f%%', colors=colors, shadow=True, startangle=0)
            PL.savefig(out_file_base + '.pie.pdf')

            PL.figure(figsize=[12,5])
            plt=PL.subplot(1,2,1)
            PL.title('Without Correction')
            p=PL.imshow(RV['beta0'][:],cmap=cm.RdBu,vmin=-0.6,vmax=+1,interpolation='None')
            PL.colorbar()
            plt.set_xticks([])
            plt.set_yticks([])
            PL.xlabel('gene'),PL.ylabel('gene')
            plt=PL.subplot(1,2,2)
            PL.title('With Correction')
            p=PL.imshow(RV['beta'][:],cmap=cm.RdBu,vmin=-0.6,vmax=+1,interpolation='None')
            PL.colorbar()
            plt.set_xticks([])
            plt.set_yticks([])
            PL.xlabel('gene'),PL.ylabel('gene')
            PL.savefig(out_file_base + '.cor.png')

            import GPy
            # Model optimization
            Ystd = Ycorr - Ycorr.mean(0)
            Ystd/=Ystd.std(0)
            ###vv = Ycorr.var(1)/(Ycorr.mean(1) ** 2)
            ###vv_idx = np.argsort(vv)[::-1]

            input_dim = 2                   # How many latent dimensions to use
            kern = GPy.kern.RBF(input_dim,ARD=True) # ARD kernel
            ###m = GPy.models.BayesianGPLVM(Ystd[vv_idx[0:500],],
            ###                             input_dim=input_dim, kernel=kern, num_inducing=40)
            m = GPy.models.BayesianGPLVM(Ystd,
                                         input_dim=input_dim, kernel=kern, num_inducing=40)
            m.optimize('scg', messages=0, max_iters=2000)
            DD = {}
            DD['PC1'] = m.X[:, 0]['mean']
            DD['PC2'] = m.X[:, 1]['mean']
            dd_out = h5py.File(out_file_base + '.PCA.h5f','w')
            dumpDictHdf5(DD,dd_out)
            dd_out.close()

            PL.figure(figsize=[6,5])
            m.kern.plot_ARD()
            PL.savefig(out_file_base + '.PCA.dim2.png')

            PL.figure(figsize=[6,5])
            # i_Gata3 = SP.where(geneID == '2625')
            i_Gata3 = SP.where(geneID == '2625')[0][0]
            PL.scatter(m.X[:,0]['mean'], m.X[:,1]['mean'], 40)
            # PL.scatter(m.X[:,0]['mean'], m.X[:,1]['mean'], 40, Y[i_Gata3, :])
            # PL.scatter(m.X[:,0]['mean'], m.X[:,1]['mean'])
            PL.xlabel('PC1')
            PL.ylabel('PC2')
            # PL.colorbar()
            PL.savefig(out_file_base + '.PCA.scatter.pdf')
