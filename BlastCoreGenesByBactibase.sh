cat only_bacteriocin_faa/* | blastp -subject BACTIBASE_renamed.faa -outfmt 6 -qcov_hsp_perc 50 | awk '{ if ( > 50) { print } }' > bacteriocins_vs_bactibase.blastp
