# enrich
Pipeline to design and test probes for target sequencing.

The first part of the pipeline is a bash script that extract the flanking sequences of target SNPs and produces 500bp fragments. Then a second bash scripts performs sliding window to produce 120 bp probes. The third part is an R script that performs BLAST and thermodynamic calculations to select probes for multiplexing.

Key packages:
rBLAST
https://github.com/mhahsler/rBLAST
primer3
https://github.com/jensenlab/primer3
