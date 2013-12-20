#!/bin/sh
#
# Output generation script for snp_meth_store

SNP_METH_STORE=../../../../../../build/Release/bin/snp_meth_store

# ============================================================
# First Section
# ============================================================


${SNP_METH_STORE} -nec -umq -mc 6 -mmq 1 -msc 5 -mpc 0.5 -if 1 -nsp -o callsNoRe_1.vcf hg18_chr21_3000.fa reads_mason2_msc5_realign_N6000.CT_GA.verified_so.sam > snp_meth_store_NoRe_1.stdout


${SNP_METH_STORE} -nse -nsi -nsd -nec -re -dr 0.01 -der 0.001 -ier 0.001 -egs 6.0 -sl -10.0 -pws 500 -mc 6 -mmq 1 -msc 5 -mpc 0.5 -if 1 -nsp -o callsRe_1.vcf hg18_chr21_3000.fa reads_mason2_msc5_realign_N6000.CT_GA.verified_so.sam  > snp_meth_store_Re_1.stdout

