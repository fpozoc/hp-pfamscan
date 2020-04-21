# High-Performance PfamScan

The objective of this pure Python implementation of PfamScan is to parallelize the process of `pfam_scan.pl` in order to perform a complete proteome.

## Installation instructions

Run the silent installation of Miniconda in case you don't have this software in your Linux Environment

```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
```

Once you have installed Miniconda/Anaconda, create a Python 3 environment.

```sh
conda create --name pfamscan python=3
conda activate pfamscan
conda install -c bioconda pfam_scan
```

In case the user does not choose Conda as the desired environment, this instructions described [here](https://vcru.wisc.edu/simonlab/bioinformatics/programs/install/pfamscan.htm) can be followed.

Disclaimer: Pfam-B has not been uploaded from version 27. You can take Pfam-A.hmm from `current_release` and Pfam-B.hmm from version 27 or take only Pfam-A.hmm. More info [here](https://en.wikipedia.org/wiki/Pfam).

```sh
git clone https://github.com/fpozoc/hp-pfamscan.git
cd HP-PfamScan
mkdir -p pfam_db

curl ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz | gunzip > pfam_db/Pfam-A.hmm
curl ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz | gunzip > pfam_db/Pfam-A.hmm.dat
curl ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz | gunzip > pfam_db/active_site.dat

curl ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-B.hmm.gz | gunzip > pfam_db/Pfam-B.hmm
curl ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-B.hmm.dat.gz | gunzip > pfam_db/Pfam-B.hmm.dat

hmmpress Pfam-A.hmm

### Optional
hmmpress Pfam-B.hmm

### Download GRCh38 Gencode v33 
mkdir -p genomes/GRCh38/g33/
curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.pc_transcripts.fa.gz | gunzip > genomes/GRCh38/g33/proteins.fa
```

Split the multifasta file in several files with transcript id as name of the file. It will be stored in `genomes/GRCh38/g33/seqs`. 

Once we have it, run `src/pfamscan.py` to locally process the sequences in a batched way.

```sh
python -m src.splitfa --infile genomes/GRCh38/g33/proteins.fa --outdir genomes/GRCh38/g33/seqs
python -m src.pfamscan --fastadir genomes/GRCh38/g33/seqs --pfamdb pfam_db
```

## Links of interest

- Some old pfam_scan.pl starting instructions [here](https://gist.github.com/olgabot/f65365842e27d2487ad3).
