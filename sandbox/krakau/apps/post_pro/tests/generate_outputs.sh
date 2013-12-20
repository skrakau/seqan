#!/bin/sh
#
# Output generation script for post_pro

POST_PRO=../../../../../../build/Release/bin/post_pro

# ============================================================
# First Section
# ============================================================

${POST_PRO} -gas -4.0 -ges -1.0 -der 0.001 -bsc 1.0 -gmr 0.2 -mq 0 -e3 4 -e4 5 -o reads_mason2_msc5_realign_N6000_simple.CT_GA.verified.sam reads_mason2_msc5_realign_N6000.CT_GA.sam hg18_chr21_3000.fa reads_mason2_msc5_realign_N6000.fastq > output1.txt


${POST_PRO} -nse -nsi -nsd -gas -4.0 -ges -1.0 -der 0.001 -bsc 1.0 -gmr 0.2 -mq 0 -e3 4 -e4 5 -o reads_mason2_msc5_realign_N6000_nonSimple.CT_GA.verified.sam reads_mason2_msc5_realign_N6000.CT_GA.sam hg18_chr21_3000.fa reads_mason2_msc5_realign_N6000.fastq > output2.txt



