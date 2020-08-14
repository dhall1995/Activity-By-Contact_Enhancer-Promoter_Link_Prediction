Due to a slight lack of datasources currently available for the cell-types/culture conditions in which I am interested, this repo will try to execute a minimally modified version of the Activity-By-Contact (ABC) model presented in https://www.nature.com/articles/s41588-019-0538-0#Sec1 for predicting enhancer promoter functional contacts based off of Hi-C and ChIPseq data. 

BACKGROUND:
During the onset of differentiation, mouse embryonic stem cells (mESCs) have been shown to pass through intermediate states of pluripotency which have different functional capacities as cells commit to particulaer lineages. The very earliest of these intermediate states is thought to be the intermediate state formed between the so-called 'naive' and 'primed' pluripotency stages. While naive stem cells are able to functionally integrate into chimeras and contribute to all lineages within the body as well as self-renew indefinitely, primed cells represent a more functionally commited group of cells with reduced ability to contribute to all lineages. Furthermore, conversion efficiency for the primed-naive transition is extremely low. Between these two states there has been shown to exist an intermediate state known as 'formative' pluripotent stem cells. These cells are thought to be the earliest non-naive mESC and are functionally different in terms of gene expression and self-renewal capacity. 

It has been previously shown by multiple sources that formative cells seem to have a very interesting state of chromatin involving changes in dynamics of methyaltion, heterochromatin formation, distribution of histone marks and Polycomb Repressive Complex (PRC) binding. Changes in the chromatin landscape are particularly evident in the distribution of different types of enhancers along the genome as well as the large growth in the number of bivalent promoters evident during the transition. Since enhancers are thought to regulate gene expression through physical interactions with promoters, it would therefore by highly useful to characterise putative enhancer-promoter interactions in naive, formative and primed mESCs to identify interactions which are gained/lost during the first stages of the exit from naive pluripotency. Such a list could be used to identify and/or test models which attempt to explain the coordinated transcriptional changes which are observed over this timecourse. 

Ideally, to identify interactions such as this we would use the previously published model (linked above) which requires both DHS and K27ac datatracks genome-wide for all the developmental timepoints in which we are interested. Since my lab (Laue Lab, Dept of Biochemistry, University of Cambridge) specialises in single-cell Hi-C a future extension would also be to incorporate these single-cell datasets for a more fine scale view of enhancer-promoter interactions. However, DHS datatracks aren't currently available for all three (naive, primed and formative) developmental timepoints and this module therefore aims to execute a minimally modified version of the ABC model using both k27ac ChIPseq signal and population Hi-C.

File formats:

- enhancer regions (e.g. enhancer_regions.bed) 
- promoter regions (e.g. promoter_regions.bed)
- gene expression (e.g. expression.bed)
- K27ac ChIPseq (BigWig format)
- Population Hi-C (.ncc format)

scripts:

- ncc_to_normalised_hic.py - takes in a .ncc format hi-c file containing the mapped reads from a hi-c experiment and creates (per-chromosome) a Hi-C map which is ~KW-normalised. We don't have incredibly high read depth in our current population Hi-C datasets so we will instead row-normalise by the square root of the sum of values in a row (following Rao et al. 2014). Optionally outputs the result to a dictionary of sparse (COO) format matrices
- filter_promoters_by_expression.py - given some expression scores and a full list of genes/promoters, returns a BED file detailing the promoters (and promoter regions + associated genes etc.) with expression higher than some threshold 
- enhancers_from_chip.py - uses ChIPseq for H3K27ac/K27me3 to identify active enhancer regions for a given timepoint
- ABC_EP_links.py - inputting ~KW-normalised per-chromosome sparse hi-C matrices, expressed promoters and enhancers (weighted by their activity score), returns the ABC links above some threshold. 