#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" src/pfamscan.py

Usage: python -m src/pfamscan.py 
--fastadir {directory with fasta splitted inside} 
--pfamdb {pfam database previously constructed}
___
--help | -h Display documentation

Example: python -m src.pfamscan --fastadir genomes/GRCh38/g27/seqs --pfamdb /media/hdd1/fpozoc/databases/Pfam/Pfam-a-33.0

Paralellized PfamScan

This standalone script aims to run pfam scan from a multisplitted fasta directory
in a batched and optimized way, using most of the processors of the machine. 

This file can also be imported as a module and contains the following functions:
    * run_pfam_scan

TO DO:  
    *
"""

from __future__ import absolute_import, division, print_function
import argparse, glob, os, sys
import functools
import multiprocessing as mp
import warnings
warnings.filterwarnings('ignore')

__author__ = "Fernando Pozo"
__copyright__ = "Copyright 2020"
__license__ = "GNU General Public License"
__version__ = "1.0.1"
__maintainer__ = "Fernando Pozo"
__email__ = "fpozoc@cnio.es"
__status__ = "Production"


def run_pfam_scan(fasta, pfamdb):
    """Function to running pfam scan from Python

    Arguments:
        fasta {str} -- Fasta file to run pfam_scan.pl
        pfamdb {str} -- Pfam database used to run the program.
    """    
    pfam_dir = f'{os.path.dirname(os.path.dirname(fasta))}/pfam_profiles'
    os.system(f'mkdir -p {pfam_dir}')
    pfam = f'{pfam_dir}/{os.path.basename(os.path.splitext(fasta)[0])}.pfam' # double iteration
    cmd = f'pfam_scan.pl -fasta {fasta} -dir {pfamdb} -outfile {pfam}'
    os.system(cmd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastadir', help='Fastas directory', type=str,
                        default='/media/hdd1/fpozoc/projects/pfam/genomes/GRCh38/g33/seqs')
    parser.add_argument('--pfamdb', help='Pfam scan database', type=str,
                        default='/media/hdd1/fpozoc/databases/Pfam/Pfam-a-33.0')
    args = parser.parse_args()

    cpus = mp.cpu_count() - 10 
    # globdir = '/media/hdd1/fpozoc/projects/pfam/genomes/GRCh38/g33/split'

    with mp.Pool(processes=cpus) as process:
        process.map(functools.partial(run_pfam_scan, pfamdb=args.pfamdb), glob.glob(f'{args.fastadir}/*'))


if __name__ == '__main__':
    main()