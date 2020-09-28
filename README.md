# Find-pBGCs
A pipeline of scripts to extract and analyze phage encoded Biosynthetic Gene Clusters (pBGCs)

# Depends on

    ncbi-genome-download: https://github.com/kblin/ncbi-genome-download
    ProphET: https://github.com/jaumlrc/ProphET
    AntiSMASH: https://github.com/antismash/antismash
    genbank_to_fasta.py: https://github.com/Coaxecva/GenBank-to-FASTA
    bioawk
    blast
    

Note that this takes very long to run, 2 days on a 64 core machine with 512gb RAM and results in ~500gb of disc use. Likely weeks, if at all, on a desktop.

# WORKFLOW

ncbi-genome-download --parallel 64 --format fasta,gff --assembly-level complete   bacteria 

Find-pBGCs.sh refseq/bacteria/
