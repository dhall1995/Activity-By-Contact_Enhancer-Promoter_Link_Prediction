##################README NOT UP TO DATE - WILL BE UPDATED SOON#############################

Due to a slight lack of datasources currently available for the cell-types/culture conditions in which I am interested, this repo will try to execute a minimally modified version of the Activity-By-Contact (ABC) model presented in https://www.nature.com/articles/s41588-019-0538-0#Sec1 for predicting enhancer promoter functional contacts based off of Hi-C and ChIPseq data. 

\subsection{Requirements}

- numpy
- pandas
- itertools
- cython
- straw

\subsection{Required input files}:

- Transcript Start Sites (TSS) (.tsv with columns ['chrom','transcrip_start','transcript_end','strand','id'] detailing the genomic positions of each TSS we are profiling) 
- H3K27ac ChIPseq file(s) (.narrowPeak format files for the approximation of regulatory element regions and assessment of their strength)

\subsection{Optional input files}: 
- Hi-C data (.hic format generated using either Juicer https://github.com/aidenlab/juicer/wiki or Pre https://github.com/aidenlab/juicer/wiki/Pre. NOTE: This algorithm assumes sufficient unique contacts to bin the genome at 5kb resolution)

\subsection{scripts}:

- make_promoter_info.py - uses the TSS inputs to create promoter regions.
- make_candidate_regions.py - Generates putative candidate elements and their strengths using the promoter_info produced using make_promoter_info.py as well as any number of H3K27ac narrowPeak files.
- linear_links.py - creates master list of possible enhancer promoter links by linking those promoters with regulatory that fall within some linear threshold distance along the chromatin backbone. Defaults to 5Mb. 
- ABC_score.py - inputting ~KW-normalised per-chromosome sparse hi-C matrices, promoter info, enhancer info, promoter expression and enhancer k27ac scores, calculates and stores a file detailing enhancer-promoter links with an ABC score above a given threshold from the master list of linear enhancer-promoter links. 

\subsection{suggested directory structure}:

/topdir/
/topdir/make_promoter_info.py
/topdir/make_candidate_regions.py
/topdir/linear_links.py
/topdir/ABC_score.py


/topdir/utils/


/topdir/data/
/topdir/data/processed/
/topdir/data/processed/links/

/topdir/data/raw/
/topdir/data/raw/input.hic
/topdir/data/raw/TSS_sites.tsv
/topdir/data/raw/tracks/{narrowPeak files}

\subsection{Example Usage:}
With the directory structure above an example usage would be:
python make_promoter_info.py --tss_sites /topdir/data/raw/TSS_sites.tsv --outpath /topdir/data/processed/

python make_candidate_elements.py --tss_sites /topdir/data/raw/TSS_sites.tsv --k27ac_path /topdir/data/raw/tracks/ --extension narrowPeak --outpath /topdir/data/processed/

python linear_links.py --regions1 /topdir/data/processed/promoter_info.csv --regions2 /topdir/data/processed/candidate_elements.csv --outpath data/processed/links/ --region1names promoter --region2names regulatory_element --distance_thresh 5000

python ABC_score.py --promoters /topdir/data/processed/promoter_info.csv --regelements /topdir/data/processed/candidate_elements.csv --k27ac /topdir/data/raw/tracks/example.narrowPeak --linearlinks /topdir/data/processed/links/linear_promoter_regulatory_element_links.npz --contacts /topdir/data/raw/input.hic --binsize 5000 --normalisation KR --threshold 0.1 --outpath /topdir/data/processed/links 
