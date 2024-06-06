This repository is linked to the article: A first passage model of intravitreal drug delivery and residence time - influence of ocular geometry, individual variability, and injection location (currently on arxiv at https://arxiv.org/abs/2404.04086).

In the folder Human_ensemble,
- The file `set_human_ensemble.ipynb` compiles the geometries used for the human ensemble, and produces Figure 4 in the manuscript.
- The file `mfpt_human_ensemble.ipynb` produces the MFPT results for the human ensemble geometries, corresponding to Figure 7b in the manuscript and Figure D.1 in the SI.
- The file `vr_exit_human_ensemble.ipynb` produces the results of the proportion of drug exiting through the vitreous-retina interface, corresponding to Figure 9 in the manuscript.

In the main manuscript, 
- Figure 7a was produced using figure_mfpt_b2_fab_igg.ipynb

In the Supplementary Information, 
- Figure A.2 was produced using figure_density_over_time.ipynb, and using the data in Data/drug_qty_over_time_P2.csv
