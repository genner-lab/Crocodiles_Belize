Crocodile hybridization paper

Data and code underpinning the paper "Population genetic structure of Morelet’s and American crocodiles in Belize: hybridization, connectivity and conservation" by Clare J. Wilkie, Marisa M. Tellez, Gareth Jones and Martin J. Genner

assets

FILE: params-crocodile.txt
The ipyrad parameter file. Note this relies on the C. porosus reference assembly GCF_001723895.fna.gz which is available at:
https://www.ncbi.nlm.nih.gov/assembly/GCF_001723895.1/

FILE: sample_codes.txt
Sample code identifiers, and the index and adaptors used for each sample.

FILE: Croc83.vcf
ipyrad output variant file, 83 specimens including reference specimens.

FILE: Croc83_80.vcf
filtered variant file, retaining only positions represented in 80% of individuals, 83 specimens including reference specimens.

FILE: Croc80_80.vcf
filtered variant file, retaining only positions represented in 80% of individuals, 80 specimens excluding 3 reference specimens.

scripts

FILE: crocodile.R
The R code used for the analyses of population genetic structure

FILE: croc_PCA_ADMIXTURE.txt
Script used for the PCA and admixture analysis
