#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" src/utils.py

Utils functions to be used in hpseqs.

This file can also be imported as a module and contains the following functions:
    * load_fasta - returns a pandas DataFrame with FASTA sequences.
    * split_multifasta  

TO DO:  
    *
"""

from __future__ import absolute_import, division, print_function

import os 
import gzip
import warnings

from Bio import SeqIO
import pandas as pd 

__author__ = "Fernando Pozo"
__copyright__ = "Copyright 2020"
__license__ = "GNU General Public License"
__version__ = "1.0.1"
__maintainer__ = "Fernando Pozo"
__email__ = "fpozoc@cnio.es"
__status__ = "Development"

warnings.filterwarnings('ignore')


def load_fasta(filepath:str)-> list:
    """Fasta Loader

    It takes a fasta file from GENCODE annotation and returns a DataFrame.
    It uses SeqIO module from BioPython (https://biopython.org/wiki/SeqIO) to parse an standard fasta
    file and retrive id and seq and parse info from the id.

    Arguments:
        filepath {str} -- fasta file path downloaded from source (.gz allowed).

    Returns:
        list -- DataFrame which contains fasta sequences, lengh and id from annotation source.
    """    
    seqlist = list()
    if filepath.endswith('.gz'):
        open_f = gzip.open(filepath, 'rt')
    else: 
        open_f = open(filepath, 'r')
    with open_f as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ids = record.name.split('|')
            seqlist.append({
                'protein_id': ids[0], 
                'transcript_id': ids[1], 
                'gene_id': ids[2],
#                 'havana_gene': ids[3],
#                 'havana_transcript': ids[4],
                'transcript_name': ids[5],
                'gene_name': ids[6],
                'id': record.name,
                'length': int(len(str(record.seq))), 
                'sequence': str(record.seq)
            })             
    df = pd.DataFrame(seqlist)
    return df


def prepare_pfamscan_outfile(filepath:str)-> list:
    """[summary]

    Arguments:
        filepath {str} -- [description]

    Returns:
        list -- [description]
    """    
    headers = ['seq_id' , 'alignment_start', 'alignment_end', 'envelope_start', 
        'envelope_end', 'hmm_acc', 'hmm_name', 'type', 'hmm_start', 
        'hmm_end', 'hmm_length', 'bit_score', 'E-value', 'significance', 'clan']
    df = pd.read_csv(filepath,  comment='#', delim_whitespace=True, names=headers)
    outpath = os.path.join(os.path.dirname(os.path.dirname(filepath)), 'spade3', os.path.basename(filepath).split('.')[0] + '.pfam')
    df.to_csv(outpath, sep='\t', index=None, headers=None)
    return df


def split_multifasta(filepath:str, outdir:str, saveby:str='transcript_id'):
    """Splitting a multifasta file in a directory

    Arguments:
        filepath {str} -- Multifasta path.
        outdir {str} -- Directory desired to split the fasta files.
        
    Keyword Arguments:
        saveby {str} -- To save fasta splitted by gene or by transcript 
    (options: {'gene_id', 'transcript_id'}) (default: {'transcript_id'})
    """    
    df_fastas = load_fasta(filepath)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if saveby == 'transcript_id':
        for i in range(0, df_fastas.shape[0]):
            pid = str(df_fastas['id'][i])
            pgene = str(df_fastas['gene_id'][i])
            ptid = str(df_fastas['transcript_id'][i])
            pseq = str(df_fastas['sequence'][i])
            with open(f"{outdir}/{ptid}.fa", "w") as fafile:
                fafile.writelines('>' + pid + '\n' + pseq + '\n')

    # this part is going to be deprecated
    elif saveby == 'gene_id':
        gid_dict = df_fastas[['gene_id', 'transcript_id']].groupby('gene_id')['transcript_id'].apply(list).to_dict()
        for gid, tlist in gid_dict.items():
            gdir = os.path.join(outdir, gid)
            if not os.path.exists(gdir):
                os.makedirs(gdir)
            for tid in tlist:
                seq = df_fastas.loc[(df_fastas['gene_id'] == gid)&(df_fastas['transcript_id'] == tid)]
                pid = seq['id'].values[0]
                pgene = seq['gene_id'].values[0]
                ptid = seq['transcript_id'].values[0]
                pseq = seq['sequence'].values[0]
                openfile = os.path.join(gdir, f'{ptid}.fa')
                with open(openfile, "w") as fastafile:
                    fastafile.writelines('>' + pid + '\n' + pseq + '\n')