# Snapper RADSeq with SISRS
Robert Literman  
Mayara Matos  

A walkthrough of an analysis to develop SNPs for differentiation *Lutjanus campechanus* and *L. purpureus*

#### Acquiring data
All data for this project was downloaded from the NCBI database. RAD-Seq data are from the paper [Genomics overrules mitochondrial DNA, siding with morphology on a controversial case of species delimitation](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2018.2924).

We also leverage an existing genome assembly for *L. campechanus* from [Development and characterization of genomic resources for a non-model marine teleost, the red snapper (Lutjanus campechanus, Lutjanidae): Construction of a high-density linkage map, anchoring of genome contigs and comparative genomic analysis](https://pubmed.ncbi.nlm.nih.gov/32348345/).

#### Step 1: Read Trimming  

All reads were trimmed using [the BBMap Suite](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/), with the following template:

```
bbduk.sh ref=adapters.fa qtrim=w ktrim=r trimq=10 maq=15 minlength=50 in=<SAMPLE>_1.fastq in2=<SAMPLE>_2.fastq out=<ID>_Trim_1.fastq.gz out2=<ID>_Trim_2.fastq.gz &> out_<ID>_Trim
```

##### Trim Results:

| SISRS_ID                | SRR        | Location | Raw_Bases     | Trim_Bases    | Percent_Surviving |
|-------------------------|------------|----------|---------------|---------------|-------------------|
| LutPur_AMA_RB1432       | SRR8647801 | AMA      | 526,021,408   | 524,055,078   | 99.6%             |
| LutPur_AMA_RB1433       | SRR8647802 | AMA      | 770,671,628   | 767,469,321   | 99.6%             |
| LutPur_AMA_RB1434       | SRR8647803 | AMA      | 533,661,648   | 529,373,060   | 99.2%             |
| LutPur_AMA_RB1436       | SRR8647804 | AMA      | 840,757,844   | 837,219,729   | 99.6%             |
| LutPur_AMA_RB1438       | SRR8647805 | AMA      | 464,558,928   | 462,632,504   | 99.6%             |
| LutPur_AMA_RB1439       | SRR8647773 | AMA      | 987,661,324   | 983,558,164   | 99.6%             |
| LutPur_AMA_RB1442       | SRR8647772 | AMA      | 517,620,756   | 515,515,316   | 99.6%             |
| LutPur_AMA_RB1444       | SRR8647771 | AMA      | 454,798,272   | 452,907,169   | 99.6%             |
| LutPur_AMA_RB1445       | SRR8647770 | AMA      | 420,242,268   | 418,608,914   | 99.6%             |
| LutPur_AMA_RB1446       | SRR8647769 | AMA      | 661,272,232   | 658,248,834   | 99.5%             |
| LutPur_AMA_RB1458       | SRR8647768 | AMA      | 705,390,232   | 702,335,336   | 99.6%             |
| LutPur_AMA_RB1459       | SRR8647767 | AMA      | 552,829,844   | 550,237,397   | 99.5%             |
| LutPur_AMA_RB1460       | SRR8647766 | AMA      | 825,226,072   | 821,079,389   | 99.5%             |
| LutPur_FOR_RB1470       | SRR8647775 | FOR      | 291,538,624   | 290,191,436   | 99.5%             |
| LutPur_FOR_RB1471       | SRR8647774 | FOR      | 649,634,368   | 646,790,826   | 99.6%             |
| LutPur_FOR_RB1473       | SRR8647857 | FOR      | 357,482,736   | 355,910,654   | 99.6%             |
| LutPur_FOR_RB1474       | SRR8647858 | FOR      | 278,134,664   | 276,805,476   | 99.5%             |
| LutPur_FOR_RB1475       | SRR8647855 | FOR      | 354,623,752   | 352,967,774   | 99.5%             |
| LutPur_FOR_RB1476       | SRR8647856 | FOR      | 468,521,464   | 466,353,392   | 99.5%             |
| LutPur_FOR_RB1477       | SRR8647853 | FOR      | 199,164,820   | 198,349,227   | 99.6%             |
| LutPur_FOR_RB1478       | SRR8647854 | FOR      | 290,538,960   | 289,353,176   | 99.6%             |
| LutPur_FOR_RB1479       | SRR8647851 | FOR      | 328,873,116   | 327,420,597   | 99.6%             |
| LutPur_FOR_RB1482       | SRR8647852 | FOR      | 317,632,744   | 316,219,035   | 99.6%             |
| LutPur_FOR_RB1484       | SRR8647860 | FOR      | 449,679,208   | 447,563,686   | 99.5%             |
| LutPur_FOR_RB1490       | SRR8647861 | FOR      | 184,981,528   | 184,003,621   | 99.5%             |
| LutPur_GUA_GUA001       | SRR8647779 | GUA      | 426,269,320   | 424,056,863   | 99.5%             |
| LutPur_GUA_GUA002       | SRR8647778 | GUA      | 1,850,542,840 | 1,839,666,299 | 99.4%             |
| LutPur_GUA_GUA003       | SRR8647785 | GUA      | 388,754,916   | 386,463,523   | 99.4%             |
| LutPur_GUA_GUA004       | SRR8647784 | GUA      | 1,045,202,376 | 1,038,667,746 | 99.4%             |
| LutPur_GUA_LP10         | SRR8647837 | GUA      | 136,079,176   | 135,311,976   | 99.4%             |
| LutPur_GUA_LP8          | SRR8647752 | GUA      | 380,908,620   | 379,044,437   | 99.5%             |
| LutPur_GUA_LP9          | SRR8647753 | GUA      | 441,330,156   | 438,975,809   | 99.5%             |
| LutPur_NUE_LP1          | SRR8647838 | NUE      | 530,394,164   | 527,690,325   | 99.5%             |
| LutPur_NUE_LP2          | SRR8647750 | NUE      | 128,357,064   | 127,571,009   | 99.4%             |
| LutPur_NUE_LP3          | SRR8647751 | NUE      | 936,876,088   | 931,885,674   | 99.5%             |
| LutPur_NUE_LP4          | SRR8647748 | NUE      | 16,810,248    | 16,378,580    | 97.4%             |
| LutPur_NUE_LP5          | SRR8647749 | NUE      | 276,921,720   | 275,434,874   | 99.5%             |
| LutPur_NUE_LP6          | SRR8647754 | NUE      | 251,656,468   | 250,385,272   | 99.5%             |
| LutPur_NUE_LP7          | SRR8647755 | NUE      | 58,902,776    | 58,562,082    | 99.4%             |
| LutPur_NUE_UMSNH16756   | SRR8647719 | NUE      | 719,661,588   | 715,661,508   | 99.4%             |
| LutPur_NUE_UMSNH16996   | SRR8647720 | NUE      | 275,952,328   | 274,214,645   | 99.4%             |
| LutPur_NUE_UMSNH16997   | SRR8647721 | NUE      | 228,268,596   | 226,999,040   | 99.4%             |
| LutPur_NUE_VEN001       | SRR8647743 | NUE      | 618,155,444   | 614,831,049   | 99.5%             |
| LutPur_NUE_VEN002       | SRR8647742 | NUE      | 772,149,280   | 768,037,575   | 99.5%             |
| LutPur_NUE_VEN003       | SRR8647741 | NUE      | 848,549,960   | 843,975,760   | 99.5%             |
| LutPur_NUE_VEN004       | SRR8647740 | NUE      | 150,922,088   | 150,146,847   | 99.5%             |
| LutPur_NUE_VEN005       | SRR8647739 | NUE      | 688,885,284   | 685,140,374   | 99.5%             |
| LutPur_NUE_VEN006       | SRR8647738 | NUE      | 323,565,368   | 321,763,036   | 99.4%             |
| LutPur_SAL_RB1520       | SRR8647828 | SAL      | 258,019,780   | 256,847,146   | 99.5%             |
| LutPur_SAL_RB1522       | SRR8647821 | SAL      | 379,174,172   | 377,507,232   | 99.6%             |
| LutPur_SAL_RB1523       | SRR8647822 | SAL      | 52,148,336    | 51,898,586    | 99.5%             |
| LutPur_SAL_RB1524       | SRR8647823 | SAL      | 78,043,796    | 77,560,358    | 99.4%             |
| LutPur_SAL_RB1525       | SRR8647824 | SAL      | 35,424,776    | 35,255,181    | 99.5%             |
| LutPur_SAL_RB1526       | SRR8647819 | SAL      | 71,993,868    | 71,703,773    | 99.6%             |
| LutPur_SAL_RB1527       | SRR8647820 | SAL      | 32,396,888    | 32,034,728    | 98.9%             |
| LutPur_SAL_RB1528       | SRR8647846 | SAL      | 295,187,088   | 293,947,204   | 99.6%             |
| LutPur_SAL_RB1531       | SRR8647845 | SAL      | 1,004,065,308 | 998,655,992   | 99.5%             |
| LutPur_SAL_RB1532       | SRR8647844 | SAL      | 1,175,408,268 | 1,168,529,824 | 99.4%             |
| LutPur_SAL_RB1533       | SRR8647843 | SAL      | 511,034,532   | 508,139,781   | 99.4%             |
| LutPur_SAL_RB1534       | SRR8647850 | SAL      | 486,619,992   | 483,905,285   | 99.4%             |
| LutPur_SAO_RB1496       | SRR8647869 | SAO      | 574,168,336   | 571,263,101   | 99.5%             |
| LutPur_SAO_RB1497       | SRR8647868 | SAO      | 730,063,460   | 725,914,957   | 99.4%             |
| LutPur_SAO_RB1498       | SRR8647871 | SAO      | 506,573,884   | 503,499,279   | 99.4%             |
| LutPur_SAO_RB1499       | SRR8647870 | SAO      | 269,630,296   | 268,373,453   | 99.5%             |
| LutPur_SAO_RB1502       | SRR8647865 | SAO      | 368,803,776   | 366,993,862   | 99.5%             |
| LutPur_SAO_RB1503       | SRR8647864 | SAO      | 320,490,180   | 318,955,962   | 99.5%             |
| LutPur_SAO_RB1504       | SRR8647867 | SAO      | 96,375,040    | 95,976,308    | 99.6%             |
| LutPur_SAO_RB1505       | SRR8647866 | SAO      | 206,464,328   | 205,483,487   | 99.5%             |
| LutPur_SAO_RB1507       | SRR8647877 | SAO      | 129,165,464   | 128,404,887   | 99.4%             |
| LutPur_SAO_RB1510       | SRR8647876 | SAO      | 271,051,188   | 269,892,762   | 99.6%             |
| LutPur_SAO_RB1514       | SRR8647825 | SAO      | 275,112,624   | 273,449,997   | 99.4%             |
| LutPur_SAO_RB1515       | SRR8647826 | SAO      | 101,413,436   | 100,793,857   | 99.4%             |
| LutPur_SAO_RB1516       | SRR8647827 | SAO      | 162,018,324   | 161,128,268   | 99.5%             |
| LutCam_ALA_ABR169       | SRR8647729 | ALA      | 88,654,820    | 88,118,211    | 99.4%             |
| LutCam_ALA_ABR170       | SRR8647728 | ALA      | 584,562,124   | 580,892,901   | 99.4%             |
| LutCam_ALA_ABR34        | SRR8647813 | ALA      | 88,385,296    | 87,775,827    | 99.3%             |
| LutCam_ALA_ABR35        | SRR8647814 | ALA      | 112,158,964   | 111,406,328   | 99.3%             |
| LutCam_ALA_ABR55        | SRR8647811 | ALA      | 91,792,100    | 91,251,476    | 99.4%             |
| LutCam_ALA_ABR59        | SRR8647812 | ALA      | 71,575,908    | 71,055,293    | 99.3%             |
| LutCam_ALA_ABR60        | SRR8647809 | ALA      | 401,296,468   | 398,794,943   | 99.4%             |
| LutCam_ALA_ABR82        | SRR8647810 | ALA      | 159,361,096   | 158,286,783   | 99.3%             |
| LutCam_ALA_ABR84        | SRR8647807 | ALA      | 422,065,124   | 419,434,794   | 99.4%             |
| LutCam_ALA_ES100614HL20 | SRR8647808 | ALA      | 128,360,848   | 127,307,236   | 99.2%             |
| LutCam_ALA_ES100614HL78 | SRR8647815 | ALA      | 37,339,996    | 37,067,105    | 99.3%             |
| LutCam_ALA_ES100614HL82 | SRR8647816 | ALA      | 205,108,968   | 203,706,898   | 99.3%             |
| LutCam_ALA_ES100614HL83 | SRR8647781 | ALA      | 265,497,480   | 263,627,237   | 99.3%             |
| LutCam_ALA_ES100614HL86 | SRR8647780 | ALA      | 109,004,484   | 108,268,066   | 99.3%             |
| LutCam_ALA_ES100614HL87 | SRR8647783 | ALA      | 233,746,108   | 232,447,573   | 99.4%             |
| LutCam_ALA_ES100614HL90 | SRR8647782 | ALA      | 505,482,200   | 502,786,829   | 99.5%             |
| LutCam_ALA_ES100614HL91 | SRR8647777 | ALA      | 299,585,300   | 297,709,154   | 99.4%             |
| LutCam_ALA_ES100614HL94 | SRR8647776 | ALA      | 25,371,204    | 25,121,484    | 99.0%             |
| LutCam_APA_LC424        | SRR8647863 | APA      | 337,265,856   | 335,123,661   | 99.4%             |
| LutCam_APA_LC425        | SRR8647806 | APA      | 292,635,468   | 290,629,547   | 99.3%             |
| LutCam_APA_LC426        | SRR8647878 | APA      | 517,524,608   | 514,249,837   | 99.4%             |
| LutCam_APA_LC427        | SRR8647879 | APA      | 448,703,796   | 445,957,626   | 99.4%             |
| LutCam_APA_LC428        | SRR8647872 | APA      | 440,274,076   | 437,524,458   | 99.4%             |
| LutCam_APA_LC429        | SRR8647873 | APA      | 239,977,840   | 238,523,345   | 99.4%             |
| LutCam_APA_LC430        | SRR8647874 | APA      | 274,393,836   | 272,819,736   | 99.4%             |
| LutCam_APA_LC431        | SRR8647875 | APA      | 246,318,964   | 244,740,885   | 99.4%             |
| LutCam_APA_LC438        | SRR8647859 | APA      | 624,653,432   | 620,553,279   | 99.3%             |
| LutCam_APA_LC447        | SRR8647862 | APA      | 22,276,408    | 21,917,176    | 98.4%             |
| LutCam_APA_LC448        | SRR8647832 | APA      | 490,870,112   | 487,845,509   | 99.4%             |
| LutCam_APA_LC449        | SRR8647831 | APA      | 470,153,572   | 467,281,512   | 99.4%             |
| LutCam_APA_LC450        | SRR8647830 | APA      | 652,293,660   | 648,188,347   | 99.4%             |
| LutCam_APA_LC451        | SRR8647829 | APA      | 608,249,964   | 604,596,019   | 99.4%             |
| LutCam_APA_LC452        | SRR8647836 | APA      | 340,179,192   | 338,110,907   | 99.4%             |
| LutCam_APA_LC453        | SRR8647835 | APA      | 515,148,944   | 511,897,088   | 99.4%             |
| LutCam_APA_LC454        | SRR8647834 | APA      | 3,175,979,312 | 3,154,577,030 | 99.3%             |
| LutCam_APA_LC455        | SRR8647833 | APA      | 3,062,776,996 | 3,043,169,418 | 99.4%             |
| LutCam_CAM_RB1555       | SRR8647849 | CAM      | 277,302,700   | 276,100,747   | 99.6%             |
| LutCam_CAM_RB1556       | SRR8647848 | CAM      | 405,206,888   | 403,215,089   | 99.5%             |
| LutCam_CAM_RB1558       | SRR8647847 | CAM      | 470,921,724   | 468,760,059   | 99.5%             |
| LutCam_CAM_RB1561       | SRR8647842 | CAM      | 407,818,364   | 405,298,023   | 99.4%             |
| LutCam_CAM_RB1562       | SRR8647841 | CAM      | 529,103,476   | 526,458,750   | 99.5%             |
| LutCam_CAM_RB1567       | SRR8647760 | CAM      | 1,064,314,328 | 1,058,734,756 | 99.5%             |
| LutCam_CAM_RB1571       | SRR8647761 | CAM      | 766,331,896   | 761,281,162   | 99.3%             |
| LutCam_CAM_UMSNH41727   | SRR8647745 | CAM      | 188,442,684   | 186,897,922   | 99.2%             |
| LutCam_CAM_UMSNH41728   | SRR8647744 | CAM      | 51,533,436    | 50,613,457    | 98.2%             |
| LutCam_FLO_309151       | SRR8647709 | FLO      | 667,967,160   | 663,626,232   | 99.4%             |
| LutCam_FLO_3091510      | SRR8647708 | FLO      | 450,500,680   | 447,409,158   | 99.3%             |
| LutCam_FLO_309153       | SRR8647711 | FLO      | 435,984,224   | 433,193,524   | 99.4%             |
| LutCam_FLO_309154       | SRR8647710 | FLO      | 552,557,052   | 548,489,811   | 99.3%             |
| LutCam_FLO_309155       | SRR8647705 | FLO      | 965,835,900   | 958,823,685   | 99.3%             |
| LutCam_FLO_309156       | SRR8647704 | FLO      | 76,794,044    | 76,011,650    | 99.0%             |
| LutCam_FLO_309157       | SRR8647707 | FLO      | 471,330,568   | 468,201,406   | 99.3%             |
| LutCam_FLO_309158       | SRR8647706 | FLO      | 813,855,496   | 808,440,539   | 99.3%             |
| LutCam_FLO_309159       | SRR8647703 | FLO      | 1,221,925,496 | 1,213,615,041 | 99.3%             |
| LutCam_FLO_3121510      | SRR8647702 | FLO      | 342,840,548   | 340,688,507   | 99.4%             |
| LutCam_FLO_3121511      | SRR8647737 | FLO      | 592,263,940   | 588,278,760   | 99.3%             |
| LutCam_FLO_3121512      | SRR8647736 | FLO      | 508,879,200   | 505,517,324   | 99.3%             |
| LutCam_FLO_3121513      | SRR8647735 | FLO      | 4,339,044     | 4,275,507     | 98.5%             |
| LutCam_FLO_3121514      | SRR8647734 | FLO      | 76,723,352    | 76,221,186    | 99.3%             |
| LutCam_FLO_3121515      | SRR8647733 | FLO      | 890,805,888   | 884,757,136   | 99.3%             |
| LutCam_FLO_3121516      | SRR8647732 | FLO      | 357,616,380   | 355,015,136   | 99.3%             |
| LutCam_FLO_3121517      | SRR8647731 | FLO      | 859,824,044   | 854,064,734   | 99.3%             |
| LutCam_FLO_312159       | SRR8647730 | FLO      | 1,083,195,112 | 1,075,464,493 | 99.3%             |
| LutCam_LOU_LSU223       | SRR8647746 | LOU      | 621,587,876   | 617,498,576   | 99.3%             |
| LutCam_LOU_LSU224       | SRR8647747 | LOU      | 362,065,160   | 359,747,911   | 99.4%             |
| LutCam_LOU_LSU226       | SRR8647840 | LOU      | 640,724,596   | 636,447,258   | 99.3%             |
| LutCam_LOU_LSU227       | SRR8647839 | LOU      | 796,398,356   | 791,318,819   | 99.4%             |
| LutCam_LOU_LSU228       | SRR8647723 | LOU      | 614,550,152   | 610,291,526   | 99.3%             |
| LutCam_LOU_LSU229       | SRR8647722 | LOU      | 854,653,724   | 848,904,401   | 99.3%             |
| LutCam_LOU_LSU230       | SRR8647725 | LOU      | 58,175,560    | 57,727,672    | 99.2%             |
| LutCam_LOU_LSU231       | SRR8647724 | LOU      | 814,744,048   | 809,563,768   | 99.4%             |
| LutCam_LOU_LSU268       | SRR8647727 | LOU      | 95,206,472    | 94,597,429    | 99.4%             |
| LutCam_LOU_LSU271       | SRR8647726 | LOU      | 163,985,660   | 162,932,529   | 99.4%             |
| LutCam_LOU_LSU272       | SRR8647818 | LOU      | 73,788,172    | 73,313,635    | 99.4%             |
| LutCam_LOU_LSU273       | SRR8647817 | LOU      | 121,848,068   | 121,030,471   | 99.3%             |
| LutCam_LOU_LSU276       | SRR8647796 | LOU      | 33,489,260    | 33,265,241    | 99.3%             |
| LutCam_LOU_LSU278       | SRR8647797 | LOU      | 643,520,112   | 639,099,400   | 99.3%             |
| LutCam_LOU_LSU279       | SRR8647798 | LOU      | 459,648,844   | 456,674,387   | 99.4%             |
| LutCam_LOU_LSU280       | SRR8647799 | LOU      | 121,950,236   | 121,207,200   | 99.4%             |
| LutCam_LOU_LSU281       | SRR8647800 | LOU      | 390,098,236   | 387,527,607   | 99.3%             |
| LutCam_TAB_RB1599       | SRR8647758 | TAB      | 901,023,548   | 895,255,813   | 99.4%             |
| LutCam_TAB_RB1600       | SRR8647759 | TAB      | 711,534,416   | 707,066,261   | 99.4%             |
| LutCam_TAB_RB1601       | SRR8647764 | TAB      | 1,063,940,916 | 1,057,284,247 | 99.4%             |
| LutCam_TAB_RB1602       | SRR8647765 | TAB      | 480,973,576   | 478,054,393   | 99.4%             |
| LutCam_TAB_RB1603       | SRR8647762 | TAB      | 728,200,184   | 723,399,931   | 99.3%             |
| LutCam_TAB_RB1604       | SRR8647763 | TAB      | 980,028,996   | 973,412,006   | 99.3%             |
| LutCam_TAB_RB1605       | SRR8647756 | TAB      | 1,278,206,820 | 1,269,750,589 | 99.3%             |
| LutCam_TUX_RB1606       | SRR8647757 | TUX      | 852,338,948   | 847,917,470   | 99.5%             |
| LutCam_TUX_RB1607       | SRR8647789 | TUX      | 386,727,896   | 384,462,105   | 99.4%             |
| LutCam_TUX_RB1608       | SRR8647788 | TUX      | 585,575,376   | 582,539,873   | 99.5%             |
| LutCam_VER_RB1609       | SRR8647791 | VER      | 130,534,412   | 129,737,277   | 99.4%             |
| LutCam_VER_RB1610       | SRR8647790 | VER      | 516,955,116   | 513,888,774   | 99.4%             |
| LutCam_VER_RB1611       | SRR8647793 | VER      | 1,161,346,924 | 1,155,095,925 | 99.5%             |
| LutCam_YUC_RB1637       | SRR8647792 | YUC      | 118,583,164   | 117,984,695   | 99.5%             |
| LutCam_YUC_RB1640       | SRR8647795 | YUC      | 355,399,988   | 353,566,500   | 99.5%             |
| LutCam_YUC_RB1641       | SRR8647794 | YUC      | 109,436,032   | 108,904,394   | 99.5%             |
| LutCam_YUC_RB1642       | SRR8647787 | YUC      | 102,940,624   | 102,435,760   | 99.5%             |
| LutCam_YUC_RB1643       | SRR8647786 | YUC      | 212,123,816   | 211,082,818   | 99.5%             |
| LutCam_YUC_RB1645       | SRR8647712 | YUC      | 172,328,004   | 171,364,513   | 99.4%             |
| LutCam_YUC_RB1646       | SRR8647713 | YUC      | 176,073,476   | 175,042,297   | 99.4%             |
| LutCam_YUC_RB1647       | SRR8647714 | YUC      | 126,940,988   | 126,275,026   | 99.5%             |
| LutCam_YUC_RB1648       | SRR8647715 | YUC      | 221,488,528   | 220,342,705   | 99.5%             |
| LutCam_YUC_RB1650       | SRR8647716 | YUC      | 80,733,360    | 80,274,791    | 99.4%             |
| LutCam_YUC_RB1651       | SRR8647717 | YUC      | 198,431,412   | 197,368,993   | 99.5%             |
| LutCam_YUC_RB1652       | SRR8647718 | YUC      | 205,525,036   | 204,449,452   | 99.5%             |

