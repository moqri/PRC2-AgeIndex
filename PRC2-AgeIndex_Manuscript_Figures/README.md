# MS Figure Generation 

### The following scripts were ran in the following order:

1.	WGBS_Mapping.sh - Mapping was ran on downloaded samples (see Supp.Table 1 in manuscript). DNMTools/Methpipe was ran to generate percentage methylation files (.meth files) used for the following analysis. Users who want the .meth files (or any other files listed in scripts) generated from online datasets can send a request to moqri@stanford.edu or dsimps93@stanford.edu. Meth files generated from fibroblast samples in our study can be downloaded from GSE253985.

2.	LMR_Generation_DNMTools.sh - LMR generation using .meth files to generate a list of low methylated regions for each given tissue type.

3.	LMR_PRC2_Binding.ipynb - Generate PRC2 binding per LMR per sample. LMR coordinates and binding saved in results folder. (t = tcell, s = epidermis, fv = passaged fibroblasts (public), c_* = various cancers, l = mouse liver.)

4.	SupplementFig2/wg2array.ipynb & Run_Clock_Comparison - Run to generate published epigenetic clock predictions of WGBS data used in the manuscript.

5.	Manuscript_Figures_PRC2-AgeIndex.ipynb - Run to generate plots and results for manuscript.

6.	PRC2_Fig6_Diffbind&heatmaps.r - Run to generate LMR heatmaps and PRC2 ChIP heatmaps in Fig 6. 

7.	SupplementFig1/PRC2_Annotation_SuppFigs.R -  Run to generate supplementary figures and analysis of LMR regions. 

