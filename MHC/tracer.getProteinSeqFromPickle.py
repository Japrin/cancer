#!/usr/bin/env python

from __future__ import print_function
import pdb
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
from operator import attrgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna

def main():
    parser = argparse.ArgumentParser(description = "get cell inform from pickle file")
    parser.add_argument('dir', metavar="<DIR>", help='directory containing subdirectories for each cell to be summarised')
    parser.add_argument('--ignore_inkt', '-i', help='ignore iNKT cells ', action="store_true")
    parser.add_argument("--sample", dest="sample_id", type=str, nargs='?', default="SAMPLE", help="sample id")
    parser.add_argument("-o", "--output_prefix", metavar="STR", type=str, dest="output_prefix",
                        default="output_prefix", help="output prefix [default: %(default)s]")
    args = parser.parse_args()

    root_dir = os.path.abspath(args.dir)

    subdirectories = os.walk(root_dir).next()[1]

    pkl_dir = "filtered_TCR_seqs"
    outdir = "{}/filtered_TCR_summary".format(root_dir)

    tracer.makeOutputDir(outdir)
    out1 = open("%s.dna_seq.fa" % args.output_prefix,"w")
    out2 = open("%s.protein_seq.fa" % args.output_prefix,"w")
    for d in subdirectories:
        cell_pkl = "{root_dir}/{d}/{pkl_dir}/{d}.pkl".format(pkl_dir=pkl_dir, d=d, root_dir=root_dir)
        if os.path.isfile(cell_pkl):
            cl = pickle.load(open(cell_pkl))
            if not cl.is_empty and not (cl.is_inkt and args.ignore_inkt):
                for locus in ['A','B','D', 'G']:
                    if cl.all_recombinants[locus] is not None:
                        for recombinant in cl.all_recombinants[locus]:
                            aaseq = Seq(str(recombinant.dna_seq), generic_dna).translate()
                            seqAnn = "productive=%s in_frame=%s stop_codon=%s cdr3=%s" % (recombinant.productive, recombinant.in_frame, recombinant.stop_codon, recombinant.cdr3)
                            print(">%s %s\n%s" % ("|".join([cl.name,locus, recombinant.contig_name, recombinant.identifier]), seqAnn, recombinant.dna_seq), file = out1)
                            print(">%s %s\n%s" % ("|".join([cl.name,locus, recombinant.contig_name, recombinant.identifier]), seqAnn, str(aaseq)), file = out2)
        print("",file = out1)
        print("",file = out2)


if __name__ == '__main__':
    main()
