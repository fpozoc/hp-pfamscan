#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" src/splitfa.py

Usage: python -m src.splitfa --infile {fastafile_path} --outdir {multifasta_dir}
___
--help | -h Display documentation

Fasta Splitted

Module which contains some functions to process fasta annotation files.
It can be used in a standalone way to split a multifasta file.

** MAIN SOURCES OF HUMAN PROTEOME ANNOTATION: **
- https://www.gencodegenes.org/human/
- https://www.uniprot.org/uniprot/
- https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml

This file can also be imported as a module and contains the following functions:
    * load_fasta: returns a pandas DataFrame with fasta annotations.
    * get_genome_ids: returns a dictionary with gene ids and nested transcript ids.
    * split_multifasta

TO DO:  
    *
"""

from __future__ import absolute_import, division, print_function
import argparse, os, sys
import warnings

import gzip
from Bio import SeqIO
import pandas as pd


__author__ = "Fernando Pozo"
__copyright__ = "Copyright 2020"
__license__ = "GNU General Public License"
__version__ = "1.0.1"
__maintainer__ = "Fernando Pozo"
__email__ = "fpozoc@cnio.es"
__status__ = "Production"

warnings.filterwarnings('ignore')


def load_fasta(filepath):
    """Fasta Loader
        
    It takes a fasta file from GENCODE annotation and returns a DataFrame.
    It uses SeqIO module from BioPython (https://biopython.org/wiki/SeqIO) to parse an standard fasta
    file and retrive id and seq and parse info from the id.


    Arguments:
        filepath {str} -- fasta file path downloaded from source


    Returns:
        df {list} -- DataFrame which contains fasta sequences, lengh and id from annotation source.
    """    
    seqlist = list()
    if filepath.endswith('.gz'):
        open_f = gzip.open(filepath, 'rt')
    else: 
        open_f = open(filepath, 'r')
    with open_f as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqlist.append({'transcript_id': record.name.split('|')[1].split('.')[0], # problems with synonyms
                            'gene_id': record.name.split('|')[2].split('.')[0],
                            'proteinLen': int(len(str(record.seq))), 
                            'proteinSeq': str(record.seq), 
                            'id': str(record.id)})             
    df = pd.DataFrame(seqlist)
    # print(f'Fasta shape = {df.shape}')            
    return df


def get_genome_ids(df_fastas):
    """Getting gene id and nested transcript id in a Python dictionary.

    Arguments:
        df_fastas {list} -- DataFrame which contains fasta sequences, lengh and id from annotation source.

    Returns:
        gdict {dict} -- Dictionary with gene id and nested transcript id
    """    
    gdict = df_fastas[['gene_id', 'transcript_id']].groupby('gene_id')['transcript_id'].apply(list).to_dict()
    return gdict


def split_multifasta(filepath, outdir, saveby='transcript_id'):
    """Split a multifasta file in a directory

    Arguments:
        filepath {str} -- Multifasta path.
        outdir {str} -- Directory desired to split the fasta files.

    Keyword Arguments:
        saveby {str} -- To save fasta splitted by gene or by transcript 
    (options: {'gene_id', 'transcript_id'}) (default: {'transcript_id'})
    """    
    df_fastas = load_fasta(filepath)
    gid_dict = get_genome_ids(df_fastas)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if saveby == 'transcript_id':
        for i in range(0, df_fastas.shape[0]):
            pid = str(df_fastas['id'][i])
            pgene = str(df_fastas['gene_id'][i])
            ptid = str(df_fastas['transcript_id'][i])
            pseq = str(df_fastas['proteinSeq'][i])
            with open(f"{outdir}/{ptid}.fa", "w") as fafile:
                fafile.writelines('>' + pid + '\n' + pseq + '\n')
                
    elif saveby == 'gene_id':
        for gid, tlist in gid_dict.items():
            dir = os.path.join(outdir, gid)
            if not os.path.exists(dir):
                os.makedirs(dir)
            for tid in tlist:
                seq = df_fastas.loc[(df_fastas['gene_id'] == gid)&(df_fastas['transcript_id'] == tid)]
                pid = seq['id'].values[0]
                pgene = seq['gene_id'].values[0]
                ptid = seq['transcript_id'].values[0]
                pseq = seq['proteinSeq'].values[0]
                with open(f"{dir}/{ptid}.fa", "w") as fafile:
                    fafile.writelines('>' + pid + '\n' + pseq + '\n')


def main():
    parser = argparse.ArgumentParser()
    # parser.add_argument("--gene_id", help="Select if you want to store your files by gene_id or not",
    #                     action='store_true', default=False)
    parser.add_argument("--infile", help="Input fasta")
    parser.add_argument("--outdir", help="Output directory")
    args = parser.parse_args()
    
    split_multifasta(args.infile, args.outdir, saveby='transcript_id')

if __name__ == '__main__':
    main()
