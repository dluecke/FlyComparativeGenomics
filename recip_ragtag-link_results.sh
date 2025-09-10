#!/bin/bash

# recip_ragtag-link_results.sh collects final output from 
# makes RT_final/ output directory and writes s links to final output assemblies

mkdir -p RT_final
cd RT_final
ln -s ../pri_haps-female/second_hap1/pri/ragtag_output-round2/ragtag.scaffold.fasta Fp_RTfinal.fa
ln -s ../pri_haps-female/second_hap1/hap1/ragtag_output-round2/ragtag.scaffold.fasta Fh1_RTfinal.fa
ln -s ../pri_haps-female/second_hap2/hap2/ragtag_output-round2/ragtag.scaffold.fasta Fh2_RTfinal.fa
ln -s ../pri_haps-male/second_hap1/pri/ragtag_output-round2/ragtag.scaffold.fasta Mp_RTfinal.fa
ln -s ../pri_haps-male/second_hap1/hap1/ragtag_output-round2/ragtag.scaffold.fasta Mh1_RTfinal.fa
ln -s ../pri_haps-male/second_hap2/hap2/ragtag_output-round2/ragtag.scaffold.fasta Mh2_RTfinal.fa
ln -s ../pri_haps-female/second_hap2/pri/ragtag_output-round2/ragtag.scaffold.fasta ALT-Fp_RTfinal.fa
ln -s ../pri_haps-male/second_hap2/pri/ragtag_output-round2/ragtag.scaffold.fasta ALT-Mp_RTfinal.fa
cd ..
