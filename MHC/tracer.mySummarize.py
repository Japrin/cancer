#!/usr/bin/env python

import pdb

import matplotlib as mpl
mpl.use('pdf')
import seaborn as sns
from matplotlib import pyplot as plt
import tracer_func as tracer
import ConfigParser
import argparse
import sys
import os
import subprocess
import pipes
import glob
import shutil
import re
from collections import defaultdict, Counter
from time import sleep
import warnings
import pickle
from prettytable import PrettyTable
import networkx as nx
from operator import attrgetter

def get_cells_component(cells, sample_id):
    cells = cells.values()
    G = nx.MultiGraph()
    #initialise all cells as nodes
    for cell in cells:
        G.add_node(cell)
    #make edges:
    for i in range(len(cells)):
        current_cell = cells[i]
        comparison_cells = cells[i+1:]
        for locus in ['A','B','D', 'G']:
            #current_identifiers = current_cell.getMainRecombinantIdentifiersForLocus(locus)
            for comparison_cell in comparison_cells:
                shared_identifiers = 0
                if current_cell.all_recombinants[locus] is not None:
                    for current_recombinant in current_cell.all_recombinants[locus]:
                        current_id_set = current_recombinant.all_poss_identifiers
                        if comparison_cell.all_recombinants[locus] is not None:
                            for comparison_recombinant in comparison_cell.all_recombinants[locus]:
                                comparison_id_set = comparison_recombinant.all_poss_identifiers
                                if len(current_id_set.intersection(comparison_id_set)) > 0:
                                    shared_identifiers += 1

                #comparison_identifiers = comparison_cell.getAllRecombinantIdentifiersForLocus(locus)
                #common_identifiers = current_identifiers.intersection(comparison_identifiers)
                if shared_identifiers > 0:
                    width = shared_identifiers * 2
                    G.add_edge(current_cell, comparison_cell, locus, penwidth=width, weight = shared_identifiers)
    deg = G.degree()
    to_remove = [n for n in deg if deg[n]==0]
    #if len(to_remove) < len(G.nodes()):
    #    G.remove_nodes_from(to_remove)
    components = nx.connected_components(G)

    out = {}
    i = 0
    print "Cell_Name\tClonotype_ID\tComponent_Size" + \
        "\tIdentifier(Alpha1)\tCDR3(Alpha1)\tTPM(Alpha1)\tProductive(Alpha1)\tIn_Frame(Alpha1)\tStop_Codon(Alpha1)" + \
        "\tIdentifier(Alpha2)\tCDR3(Alpha2)\tTPM(Alpha2)\tProductive(Alpha2)\tIn_Frame(Alpha2)\tStop_Codon(Alpha2)" + \
        "\tIdentifier(Beta1)\tCDR3(Beta1)\tTPM(Beta1)\tProductive(Beta1)\tIn_Frame(Beta1)\tStop_Codon(Beta1)" + \
        "\tIdentifier(Beta2)\tCDR3(Beta2)\tTPM(Beta2)\tProductive(Beta2)\tIn_Frame(Beta2)\tStop_Codon(Beta2)"

    for component in components:
        clonotype_id = "%s_C%04d" % (sample_id, i)
        i= i + 1
        component_size = len(component)
        for cell in component:
            output_line = "%s\t%s\t%d" % (cell.name, clonotype_id, component_size)
            for locus in ['A', 'B']:
                recombinants = cell.all_recombinants[locus]
                recombinants = sorted(recombinants, key=attrgetter('identifier'))
                if recombinants is not None:
                    for j in range(2):
                        if j < len(recombinants):
                            rec = recombinants[j]
                            if rec is not None:
                                output_line = output_line + "\t" +  "\t".join([rec.identifier, str(rec.cdr3), str(rec.TPM), str(rec.productive), str(rec.in_frame), str(rec.stop_codon)])
                            else:
                                output_line = output_line + "\t" +  "\t".join(["NA","NA","NA","NA","NA","NA"])
                        else:
                            output_line = output_line + "\t" +  "\t".join(["NA","NA","NA","NA","NA","NA"])
            print output_line
    return(out)

def main():
    parser = argparse.ArgumentParser(description = "Summarise set of cells with reconstructed TCR sequences")
    parser.add_argument('dir', metavar="<DIR>", help='directory containing subdirectories for each cell to be summarised')
    parser.add_argument('--keep_inkt', '-i', help='ignore iNKT cells when constructing networks', action="store_true")
    parser.add_argument("--sample", dest="sample_id", type=str, nargs='?', default="SAMPLE", help="sample id")
    args = parser.parse_args()

    root_dir = os.path.abspath(args.dir)

    cells = {}
    empty_cells = []
    NKT_cells = {}
    subdirectories = os.walk(root_dir).next()[1]

    pkl_dir = "filtered_TCR_seqs"
    outdir = "{}/filtered_TCR_summary".format(root_dir)

    tracer.makeOutputDir(outdir)

    for d in subdirectories:
        cell_pkl = "{root_dir}/{d}/{pkl_dir}/{d}.pkl".format(pkl_dir=pkl_dir, d=d, root_dir=root_dir)
        if os.path.isfile(cell_pkl):
            cl = pickle.load(open(cell_pkl))
            cells[d] = cl
            if cl.is_empty:
                empty_cells.append(d)
            if cl.is_inkt:
                NKT_cells[d] = (cl.is_inkt, cl.getMainRecombinantIdentifiersForLocus('B'))

    for cell_name in empty_cells:
        del cells[cell_name]
    if not args.keep_inkt:
        for cell_name in NKT_cells.keys():
            del cells[cell_name]
    #output cells' info
    get_cells_component(cells, args.sample_id)
    #plot clonotype sizes
    plt.figure()
    clonotype_sizes = tracer.get_component_groups_sizes(cells)
    w = 0.85
    x_range = range(1, len(clonotype_sizes) + 1)
    plt.bar(x_range, height=clonotype_sizes, width=w, color='black', align='center')
    plt.gca().set_xticks(x_range)
    plt.savefig("{}/clonotype_sizes_test.pdf".format(outdir))

if __name__ == '__main__':
    main()
