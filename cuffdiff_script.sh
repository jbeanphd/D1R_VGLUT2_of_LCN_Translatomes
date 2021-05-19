#! /bin/sh
cuffdiff -p 16  --no-update-check -o cuffdiff_output_folder -b Bowtie2Index/genome.fa -u genes.gtf -M exclude.gtf -L Vglut2_N,Vglut2_IP,Drd1_N,Drd1_IP Vglut2/tophat/vg1_n/accepted_hits.bam,Vglut2/tophat/vg3_n/accepted_hits.bam,Vglut2/tophat/vg4_n/accepted_hits.bam,Vglut2/tophat/vg5_n/accepted_hits.bam,Vglut2/tophat/vg6_n/accepted_hits.bam,Vglut2/tophat/vg7_n/accepted_hits.bam Vglut2/tophat/vg1_ip/accepted_hits.bam,Vglut2/tophat/vg3_ip/accepted_hits.bam,Vglut2/tophat/vg4_ip/accepted_hits.bam,Vglut2/tophat/vg5_ip/accepted_hits.bam,Vglut2/tophat/vg6_ip/accepted_hits.bam,Vglut2/tophat/vg7_ip/accepted_hits.bam D1R/tophat/D1R_1_N/accepted_hits.bam,D1R/tophat/D1R_2_N/accepted_hits.bam,D1R/tophat/D1R_3_N/accepted_hits.bam,D1R/tophat/D1R_F1_N/accepted_hits.bam,D1R/tophat/D1R_F2_N/accepted_hits.bam,D1R/tophat/D1R_F3_N/accepted_hits.bam D1R/tophat/D1R_1_IP/accepted_hits.bam,D1R/tophat/D1R_2_IP/accepted_hits.bam,D1R/tophat/D1R_3_IP/accepted_hits.bam,D1R/tophat/D1R_F1_IP/accepted_hits.bam,D1R/tophat/D1R_F2_IP/accepted_hits.bam,D1R/tophat/D1R_F3_IP/accepted_hits.bam 
