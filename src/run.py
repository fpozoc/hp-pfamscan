#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" src/run.py

Usage: python -m src.run
___
--help      | -h    Display documentation
--seqs      | -s    Protein reference file with sequences in fasta format (.gz files allowed)
--outdir    | -o    Output directory
--jobs      | -j    List of cpus used to run the programs in parallel

Example: 
python -m src.run_hpseqs  \
    --seqs /media/hdd1/fpozoc/data/genomes/GRCh38/g33/gencode.v33.pc_translations.fa.gz \
    --outdir out/GRCh38/g33 \
    --jobs 25

This file can also be imported as a module and contains the following functions:
    * multi_pfamscan

TO DO:
    *
"""

from __future__ import absolute_import, division, print_function

import argparse, glob, os
import multiprocessing as mp
import functools
import warnings

from .pfamscan import run_pfamscan
from .utils import split_multifasta

__author__ = "Fernando Pozo"
__copyright__ = "Copyright 2020"
__license__ = "GNU General Public License"
__version__ = "1.0.1"
__maintainer__ = "Fernando Pozo"
__email__ = "fpozoc@cnio.es"
__status__ = "Development"

warnings.filterwarnings('ignore')


def multi_pfamscan(fastas:str, outdir:str, pfamdb:str, cpus:int):
    """Splitting a fasta refernce file in tmp dir an running PfamScan in parallel.

    Arguments:
        multifasta {str} -- Fasta reference file (.gz allowed)
        outdir {str} -- Output directory to store the results
        pfamdb {str} -- Path where pfamdb has been stored
        cpus {int} -- Number of cpus used to run in parallel this program
    """    
    tmpdir = os.path.join(outdir, 'tmp_pfam')
    split_multifasta(fastas, tmpdir)
    with mp.Pool(processes=cpus) as process:
        process.map(functools.partial(run_pfamscan, outdir=outdir, pfamdb=), glob.glob(os.path.join(tmpdir, '*')))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seqs', '-s', 
                        help='Protein reference file with sequences in fasta format (.gz files allowed).', 
                        type=str, default='/media/hdd1/fpozoc/data/genome_annotation/GRCh38/g33/gencode.v33.pc_translations.fa.gz')
    parser.add_argument('--jobs', '-j', 
                        help='List of cpus used to run the programs in parallel.', 
                        type=int, default=25)
    parser.add_argument('--pfamdb', '-db', 
                        help='Path where pfamdb has been stored.', 
                        type=int, default=25)
    parser.add_argument('--outdir', '-o',
                        help='Output directory.', 
                        type=str, default='out/GRCh38/g33')
    args = parser.parse_args()
    
    multi_pfamscan(fastas=args.seqs, outdir=args.outdir, pfamdb=args.pfamdb, cpus=args.jobs)       
    
    os.system(f"rm -rf {os.path.join(args.outdir, 'tmp*')}")

if __name__ == '__main__':
    main()
