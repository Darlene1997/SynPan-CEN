![image](https://github.com/Darlene1997/SynPan-CEN/assets/67579562/2d233b76-00eb-4003-ba71-c10d778359fe)
# SynPan-CEN
A method to identify and compare the syntenic satellite pairs among centromeres, utilizing satellite phylogenetic information and pairwise edit distance.

### Dependencies:
- [bedtools](https://bedtools.readthedocs.io/en/latest/index.html)
- [BLASTN](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- [python-Levenshtein](https://pypi.org/project/python-Levenshtein/)
- [DAGchainer](https://vcru.wisc.edu/simonlab/bioinformatics/programs/dagchainer/dagchainer_documentation.html)

Specify your PATH to DAGCHAINER in [02_framework_ED_min_new_0601.py] and [03_addwith_window_blocks_0601.py] firstly.


### Input files
- Assembly files of reference and query genomes
- Annotation files of reference and query centromeres
```
e.g. NIP_full_annotation.rmdup.txt

Chr01_NIP	17482527	17482693	+	Chr01_NIP_M966_F4.2748	group4
Chr01_NIP	17482694	17482860	+	Chr01_NIP_M967_F4.2749	group4
Chr01_NIP	17482861	17483026	+	Chr01_NIP_M968_F4.2750	group4
Chr01_NIP	17483027	17483181	+	Chr01_NIP_M969_F6.2751	group6
Chr01_NIP	17483182	17483348	+	Chr01_NIP_M970_SF4.2752	4+
Chr01_NIP	17483349	17483514	+	Chr01_NIP_M971_F4.2753	group4
Chr01_NIP	17483515	17483616	+	Chr01_NIP_else_SF14.2754	14+
Chr01_NIP	17483617	17483723	+	Chr01_NIP_else_SF2.2755	2+
Chr01_NIP	17483724	17483887	+	Chr01_NIP_M972_F4.2756	group4
Chr01_NIP	17484123	17485044	-	Gypsy-215_OS-LTR	LTR
Chr01_NIP	17486458	17487251	-	SZ-22_LTR	LTR
Chr01_NIP	17488157	17489008	+	LOC_Os01g30970.1.MSUv7.0	Gene
Chr01_NIP	17489009	17490854	-	CW06_rnd-5_family-1612	LTR
Chr01_NIP	17490883	17490952	+	CW06_rnd-5_family-7862	other_repeat
Chr01_NIP	17491002	17491068	-	NH285_rnd-4_family-4525	other_repeat
```
- Sequences of all satellite repeats
```
e.g. CR_CENtype_repeats_NIP.csv

16810093.0,16810247.0,155,CATGATTTTTGGACATTTTGGATTGTATTGGGTGCGTTCGTGGCAAAAACTCACTTCGTGATTCGCGCGGCGAACTTTTGTTTATTAATGCGAATATTGGCACACGAGGGTGCGATGTTTTTGACCATAATCAAAAAGTTTAAAAAACGCCAAAA,-,CEN154,NIP.fasta_Chr01_NIP,Chr01_NIP,13.0,nan
16810248.0,16810413.0,166,CATGATTTTTGGACATATTGGAGTGTATTGGGTGTGTTCGTGGCAAAAACTCACTTCGTGAATTGCGCGACGAAATTTTGTCAATTAATGCCAATATTGCTATATTTTGGCACACACGGGTGCGTTTTTTTTTGTACCGGAATGAAAAGTTCGAAAAGCACCAAAA,-,CEN154,NIP.fasta_Chr01_NIP,Chr01_NIP,21.0,nan
16810414.0,16810578.0,165,CATGATTTTTGGACATATTGGAGTGTATTGGGTGTGTTCGTGGCAAAAACTCACTTCGTGATTCGTGCGGCGAACTCTTGTCAAATAATGCCAATATTGGCATATTTTGGCCCGCACAGGTGCGATGTTTTTGACCGGAATGAAAAAGTTCAAAAAGCATCAAAA,-,CEN154,NIP.fasta_Chr01_NIP,Chr01_NIP,17.0,nan
16810579.0,16810743.0,165,CATGATTTTTGGACATATTGGAGTGTATTGGGTGCGTTCATGGCAAAAACTCACTTCGTGATTCGCGCGGCGAAATTTTGTCAATTAATTCCAATGTTGGCATATTTTGGCCCACACGGTTGCGATGTTTTTGACCGGAATGAAAAAGTTCAGAAAGCGACAAAA,-,CEN154,NIP.fasta_Chr01_NIP,Chr01_NIP,18.0,nan
16810744.0,16810908.0,165,CATGATTTTTGGACATATTGGAGTGTATTGGATGTGTTCGTGGCAAAAACTCACTTGGTGATTCGCGCGGCGAACTTTTGTCACTTAATGCCAATATTGGCATATTTTCGTCCACACGGGTGCGATGTTTTTGACCGGAATGAAAAAGTTCAAAAAGCACCAAAA,-,CEN154,NIP.fasta_Chr01_NIP,Chr01_NIP,14.0,nan
16810909.0,16811073.0,165,CATGATTTTTGGACATATTGGAGTGTATTGGATGTGTTCGTGGCAAAAACTCACTTGGTGATTCGCGCGGTGTACTTTTGTCAGTTAATGCCAATATTGGCATATTTTGGCCCGTACGGGTGCGATGTTTTTAACTGGAGTGAAAAAGTTCAAAAAGCACCAAGA,-,CEN154,NIP.fasta_Chr01_NIP,Chr01_NIP,21.0,nan
16811074.0,16811238.0,165,CATGATTTTTGGACATATTGGAGTGTATTGGGTGTGTTCGTGGCAAAAACTCACTTCGTGATTCGCGCGGCGTACTTTTGTCAGTTAATGCCAATATAGGCATATTCTGGCCCGCACGGGTGCAATGTTTTTGACTGAAATGAAAAAGTTCAAAAAGCACCAAGA,-,CEN154,NIP.fasta_Chr01_NIP,Chr01_NIP,19.0,na
```
### Run
- Get blast best hit of PAIRWISE LTRs in two accessions
```
python LTR_BlastBestHit_Pairs.py -r ${ref} -q ${que} -chr Chr${each_chr} -d ${Dir}
```

- Get all-vs-all pairwise ED comparisons within the same group
```
python Monomer_PairED_InGroup.py -r ${ref} -q ${que} -chr Chr${each_chr} -d ${Dir}
```

- Find a framework of minimum edit distance
```
python Framework_WithMinED.py -r ${ref} -q ${que} -chr Chr${each_chr} -d ${Dir}
```

- Slide windows for DAGchainer to decide the local best alignments
```
python SlideWindow_DAGchainer.py -r ${ref} -q ${que} -chr Chr${each_chr} -e 15 -d ${Dir}
```

- Get the formats for ploting
```
python Merge_Tracks_Pairs.py -r ${ref} -q ${que} -c Chr${each_chr} -e 15 -d ${Dir}
python ED_ColorRange.py -r ${ref} -q ${que} -c Chr${each_chr} -e 15 -d ${Dir}
python Add_More_Details.py -r ${ref} -q ${que} -c Chr${each_chr} -e 15 -d ${Dir}
```

- Ploting
```
Rscript MatchList_Pairs_12chr_EDColorRange.R ${ref}_$(que) 15
```


### Output files
(e.g. using Nippobare as reference and SL177 as query)
- Chr01_NIP_SL177_LRT_blast_best_hit.txt
```
Chr01_NIP       Chr01_NIP.LTR.28        16814766        16815559        Chr01_SL177     Chr01_SL177.LTR.28      17288185        17288978        0.0     0
Chr01_NIP       Chr01_NIP.LTR.2443      17418410        17419203        Chr01_SL177     Chr01_SL177.LTR.28      17288185        17288978        0.0     0
Chr01_NIP       Chr01_NIP.LTR.2572      17449648        17450439        Chr01_SL177     Chr01_SL177.LTR.28      17288185        17288978        0.0     0
Chr01_NIP       Chr01_NIP.LTR.2758      17486458        17487251        Chr01_SL177     Chr01_SL177.LTR.28      17288185        17288978        0.0     0
Chr01_NIP       Chr01_NIP.LTR.785       16970834        16971626        Chr01_SL177     Chr01_SL177.LTR.28      17288185        17288978        0.0     0
Chr01_NIP       Chr01_NIP.LTR.720       16956863        16957655        Chr01_SL177     Chr01_SL177.LTR.28      17288185        17288978        0.0     0
Chr01_NIP       Chr01_NIP.LTR.717       16952045        16952836        Chr01_SL177     Chr01_SL177.LTR.28      17288185        17288978        0.0     0
```


- Chr01_NIP_SL177_group_pairED_all.matchList
```
Chr01_NIP       Chr01_NIP_M1022_F1.2814 17509781        17509945        Chr01_SL177     Chr01_SL177_M1022_F1.2965       18023283        18023447        0.0     0
Chr01_NIP       Chr01_NIP_M6_F2.4       16810744        16810908        Chr01_SL177     Chr01_SL177_M6_F2.5     17284166        17284330        0.0     0
Chr01_NIP       Chr01_NIP_M6_F2.4       16810744        16810908        Chr01_SL177     Chr01_SL177_M30_F2.32   17289438        17289602        8.484848484848485e-06   14
Chr01_NIP       Chr01_NIP_M6_F2.4       16810744        16810908        Chr01_SL177     Chr01_SL177_M32_F2.34   17289768        17289932        7.878787878787878e-06   13
Chr01_NIP       Chr01_NIP_M6_F2.4       16810744        16810908        Chr01_SL177     Chr01_SL177_M39_F2.42   17291037        17291201        7.878787878787878e-06   13
Chr01_NIP       Chr01_NIP_M6_F2.4       16810744        16810908        Chr01_SL177     Chr01_SL177_M153_F2.127 17316050        17316214        6.666666666666667e-06   11
Chr01_NIP       Chr01_NIP_M6_F2.4       16810744        16810908        Chr01_SL177     Chr01_SL177_M159_F2.134 17317132        17317296        6.666666666666667e-06   11
Chr01_NIP       Chr01_NIP_M6_F2.4       16810744        16810908        Chr01_SL177     Chr01_SL177_M978_F2.2918        18015047        18015211        7.878787878787878e-06   13
Chr01_NIP       Chr01_NIP_M153_F2.386   16882118        16882282        Chr01_SL177     Chr01_SL177_M6_F2.5     17284166        17284330        6.666666666666667e-06   11
Chr01_NIP       Chr01_NIP_M153_F2.386   16882118        16882282        Chr01_SL177     Chr01_SL177_M30_F2.32   17289438        17289602        7.878787878787878e-06   13
Chr01_NIP       Chr01_NIP_M153_F2.386   16882118        16882282        Chr01_SL177     Chr01_SL177_M32_F2.34   17289768        17289932        7.272727272727272e-06   12
Chr01_NIP       Chr01_NIP_M153_F2.386   16882118        16882282        Chr01_SL177     Chr01_SL177_M39_F2.42   17291037        17291201        7.272727272727272e-06   12
Chr01_NIP       Chr01_NIP_M153_F2.386   16882118        16882282        Chr01_SL177     Chr01_SL177_M153_F2.127 17316050        17316214        0.0     0
Chr01_NIP       Chr01_NIP_M153_F2.386   16882118        16882282        Chr01_SL177     Chr01_SL177_M159_F2.134 17317132        17317296        6.060606060606061e-0

```


- Chr01_NIP_SL177_group_pairED_order_all.dag (This file is used for DAGchainer)
```
Chr01_NIP       Chr01_NIP_M1022_F1.2814 2814    2814    Chr01_SL177     Chr01_SL177_M1022_F1.2965       2965    2965    1e-50   0
Chr01_NIP       Chr01_NIP_M6_F2.4       4       4       Chr01_SL177     Chr01_SL177_M6_F2.5     5       5       1e-50   0
Chr01_NIP       Chr01_NIP_M6_F2.4       4       4       Chr01_SL177     Chr01_SL177_M30_F2.32   32      32      1e-20   14
Chr01_NIP       Chr01_NIP_M6_F2.4       4       4       Chr01_SL177     Chr01_SL177_M32_F2.34   34      34      1e-20   13
Chr01_NIP       Chr01_NIP_M6_F2.4       4       4       Chr01_SL177     Chr01_SL177_M39_F2.42   42      42      1e-20   13
Chr01_NIP       Chr01_NIP_M6_F2.4       4       4       Chr01_SL177     Chr01_SL177_M153_F2.127 127     127     1e-20   11
Chr01_NIP       Chr01_NIP_M6_F2.4       4       4       Chr01_SL177     Chr01_SL177_M159_F2.134 134     134     1e-20   11
Chr01_NIP       Chr01_NIP_M6_F2.4       4       4       Chr01_SL177     Chr01_SL177_M978_F2.2918        2918    2918    1e-20   13
Chr01_NIP       Chr01_NIP_M153_F2.386   386     386     Chr01_SL177     Chr01_SL177_M6_F2.5     5       5       1e-20   11
Chr01_NIP       Chr01_NIP_M153_F2.386   386     386     Chr01_SL177     Chr01_SL177_M30_F2.32   32      32      1e-20   13
Chr01_NIP       Chr01_NIP_M153_F2.386   386     386     Chr01_SL177     Chr01_SL177_M32_F2.34   34      34      1e-20   12
Chr01_NIP       Chr01_NIP_M153_F2.386   386     386     Chr01_SL177     Chr01_SL177_M39_F2.42   42      42      1e-20   12
Chr01_NIP       Chr01_NIP_M153_F2.386   386     386     Chr01_SL177     Chr01_SL177_M153_F2.127 127     127     1e-50   0
Chr01_NIP       Chr01_NIP_M153_F2.386   386     386     Chr01_SL177     Chr01_SL177_M159_F2.134 134     134     1e-40   1
Chr01_NIP       Chr01_NIP_M153_F2.386   386     386     Chr01_SL177     Chr01_SL177_M978_F2.2918        2918    2918    1e-30   4
```


- Chr01_NIP_SL177_group_pairED0_framework.sg
```
The name of the file means use ED=0 as the edit distance cutoff for the framework construction.

Chr01_NIP       Chr01_NIP_M2_F8.0       0       0       Chr01_SL177     Chr01_SL177_M2_F8.1     1       1       9.999999999999999e-51   50
Chr01_NIP       Chr01_NIP_M4_SF2.2      2       2       Chr01_SL177     Chr01_SL177_M4_SF2.3    3       3       9.999999999999999e-51   97
Chr01_NIP       Chr01_NIP_M5_SF2.3      3       3       Chr01_SL177     Chr01_SL177_M5_SF2.4    4       4       9.999999999999999e-51   147
Chr01_NIP       Chr01_NIP_M6_F2.4       4       4       Chr01_SL177     Chr01_SL177_M6_F2.5     5       5       9.999999999999999e-51   197
Chr01_NIP       Chr01_NIP_M7_SF5.5      5       5       Chr01_SL177     Chr01_SL177_M7_SF5.6    6       6       9.999999999999999e-51   247
Chr01_NIP       Chr01_NIP_M8_SF5.6      6       6       Chr01_SL177     Chr01_SL177_M8_SF5.7    7       7       9.999999999999999e-51   297
Chr01_NIP       Chr01_NIP_M9_SF5.7      7       7       Chr01_SL177     Chr01_SL177_M9_SF5.8    8       8       9.999999999999999e-51   347
Chr01_NIP       Chr01_NIP_M10_SF1.8     8       8       Chr01_SL177     Chr01_SL177_M10_SF1.9   9       9       9.999999999999999e-51   397
Chr01_NIP       Chr01_NIP_M11_SF1.9     9       9       Chr01_SL177     Chr01_SL177_M11_SF1.10  10      10      9.999999999999999e-51   447
Chr01_NIP       Chr01_NIP_M12_F8.10     10      10      Chr01_SL177     Chr01_SL177_M12_F8.11   11      11      9.999999999999999e-51   497
......
```


- Chr01_NIP_SL177_pairED15.sg_filter(pairwise alignment result)
```
The name of the file means use ED=15 as the edit distance cutoff for the padding synteny.

Chr01_NIP       Chr01_NIP_M2_F8.0       0       0       Chr01_SL177     Chr01_SL177_M2_F8.1     1       1
Chr01_NIP       Chr01_NIP_M9313_SF2.1   1       1       Chr01_SL177     Chr01_SL177_M213098_SF2.2       2       2
Chr01_NIP       Chr01_NIP_M4_SF2.2      2       2       Chr01_SL177     Chr01_SL177_M4_SF2.3    3       3
Chr01_NIP       Chr01_NIP_M5_SF2.3      3       3       Chr01_SL177     Chr01_SL177_M5_SF2.4    4       4
Chr01_NIP       Chr01_NIP_M6_F2.4       4       4       Chr01_SL177     Chr01_SL177_M6_F2.5     5       5
Chr01_NIP       Chr01_NIP_M7_SF5.5      5       5       Chr01_SL177     Chr01_SL177_M7_SF5.6    6       6
Chr01_NIP       Chr01_NIP_M8_SF5.6      6       6       Chr01_SL177     Chr01_SL177_M8_SF5.7    7       7
Chr01_NIP       Chr01_NIP_M9_SF5.7      7       7       Chr01_SL177     Chr01_SL177_M9_SF5.8    8       8
Chr01_NIP       Chr01_NIP_M10_SF1.8     8       8       Chr01_SL177     Chr01_SL177_M10_SF1.9   9       9
Chr01_NIP       Chr01_NIP_M11_SF1.9     9       9       Chr01_SL177     Chr01_SL177_M11_SF1.10  10      10
Chr01_NIP       Chr01_NIP_M12_F8.10     10      10      Chr01_SL177     Chr01_SL177_M12_F8.11   11      11
```


- NIP_SL177_12chr_ColorRange_ED15.txt
- NIP_SL177_12chr_for_plot_ED15.txt

These two files were used for 06_matchList_of_pairs_ingroup_chr_ED_0522.R.



- Alignment plots
```
The color range between tracks represent edit distance of the element pair.
The three tracks in each side represent monomer subfamilies, strand orientations, and TEs, respectively.
```

![image](https://github.com/Darlene1997/SynPan-CEN/assets/67579562/37d1e7c5-2794-49aa-9d5e-e33bf0f5de81)


### Others

- Calculate and compare substitution rates between centromere regions and chromosome arm
```
python ../01_ED_cal_for_chromosome.py -r NIP -q ${que} -c Chr${each_chr}
python ../02_merge_for_plot.py -r NIP -q ${que} -c Chr${each_chr}
```

If any questions, contact at xielj1122@foxmail.com.

Cite: 
