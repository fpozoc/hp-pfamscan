#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" src/pfamscan.py

This file can also be imported as a module and contains the following functions:
    * pfamscan
    * run_pfam_scan

TO DO:  
    *
"""
from __future__ import absolute_import, division, print_function

import os
import warnings

from .utils import split_multifasta

__author__ = "Fernando Pozo"
__copyright__ = "Copyright 2020"
__license__ = "GNU General Public License"
__version__ = "1.0.1"
__maintainer__ = "Fernando Pozo"
__email__ = "fpozoc@cnio.es"
__status__ = "Development"

warnings.filterwarnings('ignore')


def pfamscan(infile:str, outfile:str, pfamdb:str)-> str:
    """PfamScan
    PfamScan is used to search a FASTA sequence against a library of Pfam HMM.
    https://www.ebi.ac.uk/Tools/pfa/pfamscan/

    Install the command line interface from:
    https://anaconda.org/bioconda/pfam_scan

    Arguments:
        infile {str} -- Specified using the -fasta option. Must be in FASTA format.
        outfile {str} -- Specified using the -outfile option.
        pfamdb {str} -- directory location of Pfam files.

    Returns:
        str -- Command line argument selected.
    """    
    cmd = f'pfam_scan.pl -fasta {infile} -dir {pfamdb} -outfile {outfile}'
    return cmd 


def run_pfamscan(fasta:str, outdir:str, pfamdb:str):
    """Running `pfam_scan.pl` for a common FASTA file and storing it in pfamscan
    directory.

    Arguments:
        fasta {str} -- FASTA file to run pfam_scan.pl.
        outdir{str} -- Directory to store the results.
        pfamdb {str} -- Pfam database used to run the program.
    """    
    outdir = os.path.join(os.path.dirname(os.path.dirname(fasta)), 'pfamscan')
    os.system(f'mkdir -p {outdir}')
    outfile = os.path.join(outdir, os.path.basename(os.path.splitext(fasta)[0]) + '.pfam')    
    # pfamdb = '/media/hdd1/fpozoc/databases/Pfam/Pfam-a-33.0' # changeable
    os.system(pfamscan(fasta, outfile, pfamdb))