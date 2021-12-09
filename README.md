# Snapper RADSeq with SISRS
Robert Literman  
Mayara Matos  

A walkthrough of an analysis to develop SNPs for differentiation *Lutjanus campechanus* and *L. purpureus*

#### Acquiring data
All data for this project was downloaded from the NCBI database. Most data are from the paper [Genomics overrules mitochondrial DNA, siding with morphology on a controversial case of species delimitation](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2018.2924), but some also come from [dDocent: a RADseq, variant-calling pipeline designed for population genomics of non-model organisms](https://peerj.com/articles/431/).

#### Step 1: Assemble a Composite Genome with SISRS

To assemble the SISRS composite genome, we downloaded data from the far reaches of the range to avoid introgression issues. 18 *L. campechanus* samples came from Florida, USA, while the 12 *L. purpureus* samples came from Fortaleza, Brazil.

```
fasterq-dump SRR8647703
fasterq-dump SRR8647730
fasterq-dump SRR8647705
fasterq-dump SRR8647733
fasterq-dump SRR8647731
fasterq-dump SRR8647706
fasterq-dump SRR8647708
fasterq-dump SRR8647737
fasterq-dump SRR8647710
fasterq-dump SRR8647702
fasterq-dump SRR8647736
fasterq-dump SRR8647707
fasterq-dump SRR8647709
fasterq-dump SRR8647711
fasterq-dump SRR8647732
fasterq-dump SRR8647704
fasterq-dump SRR8647734
fasterq-dump SRR8647735
fasterq-dump SRR8647774
fasterq-dump SRR8647856
fasterq-dump SRR8647860
fasterq-dump SRR8647857
fasterq-dump SRR8647855
fasterq-dump SRR8647851
fasterq-dump SRR8647852
fasterq-dump SRR8647775
fasterq-dump SRR8647854
fasterq-dump SRR8647858
fasterq-dump SRR8647861
fasterq-dump SRR8647853
```

All reads were trimmed using [the BBMap Suite](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/), with the following template:

```
bbduk.sh ref=adapters.fa qtrim=w ktrim=r trimq=10 maq=15 minlength=50 in=<SAMPLE>_1.fastq in2=<SAMPLE>_2.fastq out=<ID>_Trim_1.fastq.gz out2=<ID>_Trim_2.fastq.gz &> out_<ID>_Trim

```

##### Trim Results for Composite Genome Samples:

| SISRS_ID           | Species | SRR        | Location | Raw_Bases     | Trim_Bases    | Percent_Surviving |
|--------------------|---------|------------|----------|---------------|---------------|-------------------|
| LutCam_FLO_309151  | LutCam  | SRR8647709 | FLO      | 667,967,160   | 663,626,232   | 99.4%             |
| LutCam_FLO_3091510 | LutCam  | SRR8647708 | FLO      | 450,500,680   | 447,409,158   | 99.3%             |
| LutCam_FLO_309153  | LutCam  | SRR8647711 | FLO      | 435,984,224   | 433,193,524   | 99.4%             |
| LutCam_FLO_309154  | LutCam  | SRR8647710 | FLO      | 552,557,052   | 548,489,811   | 99.3%             |
| LutCam_FLO_309155  | LutCam  | SRR8647705 | FLO      | 965,835,900   | 958,823,685   | 99.3%             |
| LutCam_FLO_309156  | LutCam  | SRR8647704 | FLO      | 76,794,044    | 76,011,650    | 99.0%             |
| LutCam_FLO_309157  | LutCam  | SRR8647707 | FLO      | 471,330,568   | 468,201,406   | 99.3%             |
| LutCam_FLO_309158  | LutCam  | SRR8647706 | FLO      | 813,855,496   | 808,440,539   | 99.3%             |
| LutCam_FLO_309159  | LutCam  | SRR8647703 | FLO      | 1,221,925,496 | 1,213,615,041 | 99.3%             |
| LutCam_FLO_3121510 | LutCam  | SRR8647702 | FLO      | 342,840,548   | 340,688,507   | 99.4%             |
| LutCam_FLO_3121511 | LutCam  | SRR8647737 | FLO      | 592,263,940   | 588,278,760   | 99.3%             |
| LutCam_FLO_3121512 | LutCam  | SRR8647736 | FLO      | 508,879,200   | 505,517,324   | 99.3%             |
| LutCam_FLO_3121513 | LutCam  | SRR8647735 | FLO      | 4,339,044     | 4,275,507     | 98.5%             |
| LutCam_FLO_3121514 | LutCam  | SRR8647734 | FLO      | 76,723,352    | 76,221,186    | 99.3%             |
| LutCam_FLO_3121515 | LutCam  | SRR8647733 | FLO      | 890,805,888   | 884,757,136   | 99.3%             |
| LutCam_FLO_3121516 | LutCam  | SRR8647732 | FLO      | 357,616,380   | 355,015,136   | 99.3%             |
| LutCam_FLO_3121517 | LutCam  | SRR8647731 | FLO      | 859,824,044   | 854,064,734   | 99.3%             |
| LutCam_FLO_312159  | LutCam  | SRR8647730 | FLO      | 1,083,195,112 | 1,075,464,493 | 99.3%             |
| LutPur_FOR_RB1470  | LutPur  | SRR8647775 | FOR      | 291,538,624   | 290,191,436   | 99.5%             |
| LutPur_FOR_RB1471  | LutPur  | SRR8647774 | FOR      | 649,634,368   | 646,790,826   | 99.6%             |
| LutPur_FOR_RB1473  | LutPur  | SRR8647857 | FOR      | 357,482,736   | 355,910,654   | 99.6%             |
| LutPur_FOR_RB1474  | LutPur  | SRR8647858 | FOR      | 278,134,664   | 276,805,476   | 99.5%             |
| LutPur_FOR_RB1475  | LutPur  | SRR8647855 | FOR      | 354,623,752   | 352,967,774   | 99.5%             |
| LutPur_FOR_RB1476  | LutPur  | SRR8647856 | FOR      | 468,521,464   | 466,353,392   | 99.5%             |
| LutPur_FOR_RB1477  | LutPur  | SRR8647853 | FOR      | 199,164,820   | 198,349,227   | 99.6%             |
| LutPur_FOR_RB1478  | LutPur  | SRR8647854 | FOR      | 290,538,960   | 289,353,176   | 99.6%             |
| LutPur_FOR_RB1479  | LutPur  | SRR8647851 | FOR      | 328,873,116   | 327,420,597   | 99.6%             |
| LutPur_FOR_RB1482  | LutPur  | SRR8647852 | FOR      | 317,632,744   | 316,219,035   | 99.6%             |
| LutPur_FOR_RB1484  | LutPur  | SRR8647860 | FOR      | 449,679,208   | 447,563,686   | 99.5%             |
| LutPur_FOR_RB1490  | LutPur  | SRR8647861 | FOR      | 184,981,528   | 184,003,621   | 99.5%             |

