# enRich
Pipeline to design and test probes for target sequencing.

Auxiliary scripts to select targets for enRich based on SNPeff annotations. The first part of the pipeline is a bash script that extracts the flanking sequences of target SNPs. Then a second bash script is called to performs sliding window to produce 120 bp probes. Finally BLAST is performed on the candidate probes to identify the potential binding sites and off targets. The third part is an R script that takes the BLAST results and the sequences themselves to perform filtering for uniqueness, GC, Tm and thermodynamic calculations such as haipin and homodimer to select probes for multiplexing.

Other bash script for fasta sliding window, please refer to:
https://github.com/kdillmcfarland/sliding_windows

A key component of the pipeline is BLAST. It is important to make sure that BLAST is correctly installed and added to the PATH. For more information on BLAST configuration please refer to:
https://www.ncbi.nlm.nih.gov/books/NBK569861/

Key package for R:

primer3
https://github.com/jensenlab/primer3

Biostrings
https://www.bioconductor.org/packages/release/bioc/html/Biostrings.html

data.table
https://cran.r-project.org/web/packages/data.table/index.html

dplyr
https://dplyr.tidyverse.org/

stringr
https://stringr.tidyverse.org/index.html

ggplot2
https://ggplot2.tidyverse.org/

![alt text](https://github.com/renanso/enrich_v1/blob/main/scheme.png "Logo Title Text 1")

