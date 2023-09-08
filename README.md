[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7604944.svg)](https://doi.org/10.5281/zenodo.7604944)

**Crocodile hybridization paper**

Data and code underpinning the paper "Population genetic structure of Morelet’s and American crocodiles in Belize: hybridization, connectivity and conservation" by Clare J. Wilkie, Marisa M. Tellez, Gareth Jones and Martin J. Genner

assets

FILE: params-crocodile.txt
The ipyrad parameter file. Note this relies on the C. porosus reference assembly GCF_001723895.fna.gz which is available at:
https://www.ncbi.nlm.nih.gov/assembly/GCF_001723895.1/

FILE: Croc83.vcf
ipyrad output variant file, 83 specimens including reference specimens.

FILE: Croc83_80.vcf
filtered variant file, retaining only positions represented in 80% of individuals, 83 specimens including reference specimens.

FILE: Croc80_80.vcf
filtered variant file, retaining only positions represented in 80% of individuals, 80 specimens excluding 3 reference specimens.

FILE: Admixture_Long.txt
used for plotting admixture values. sample = sample code; species and prob combine to give the proportional genetic composition of the individual to species

FILE: PCA_plot_data.txt
used for PCA plotting. sample = sample code; PC1 = values on PC1; PC2 = values on PC2; Group = group category based on admixture results for wild samples

FILE: Population_80.txt used for assigning wild individuals to populations in population genetic analysis. Species = assigned using admixture results. Population = assigned using geographical data.

scripts

FILE: CrocodileScript.R
The R code used for the analyses of population genetic structure

FILE: croc_PCA_ADMIXTURE.txt
Script used for the PCA and admixture analysis

FILE: FineRADstructure_script.txt
Script used for FineRADstructure analysis
