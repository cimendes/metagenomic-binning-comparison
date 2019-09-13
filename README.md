#  Benchmarking of taxonomic independent metagenomic binning software for short-read data

This is an atempt to benchmark with relevant metagenomic binning software that does not rely on a reference database.

:warning: WORK IN PROGRESS :warning:


## Table of Contents

* [Introduction](#introduction)
* [Methods](#methods)
    * [Metagenomic Dataset](#metagenomic-datasets)
        * [simulated dataset](#simulated-dataset)
        * [Zymos community standard](#zymos-community-standard)
        * [Real dataset](#real-dataset)
    * [Binning tools and commands](#binning-tools-and-commands)
    * [Assessing Binning Success](#assessing-binning-success)
* [Results](#resuls)
* [References](#references)
* [Authors](#authors)


## Introduction
Shotgun metagenomics can offer relatively unbiased pathogen detection and characterization, potentially able to provide genotyping, antimicrobial resistance and virulence profiling in a single methodological step. This comes with the cost of producing massive amounts of information that require expert handling and processing, as well as capable computational infrastructures. One of the biggest challenges when dealing with metagenomic data is the lack of gold standards, although major efforts are being made on the standardization and assessment of software, both commercial and open source (Angers-Loustau et al., 2018; Gruening et al., 2018; Sczyrba et al., 2017: Couto et al., 2018). 

A plethora of tools are available specifically for metagenomic data, both short and long read data, and several combinations of these tools can be used to characterize the causative agent in a patient's infection in a fraction of the time required by traditional methods. The basic strategy for analysing metagenomic data is either a direct classification and characterization of the sequencing data, a metagenomic assembly followed by taxonomic independent or dependent binning that is then classified and characterized, or a combination of both. The assembly methods provide longer sequences are more informative than shorter sequencing data and can provide a more complete picture of the microbial community in a given sample. 

In metagenomics, the grouping of the different contigs obtained in the metagenomic assembly, ideally each collecting the sequences that belong to a microorganism present in the sample, represent the greatest bottleneck when trying to obtain fiable, reproducible results. The euristics implemented in this process, necessary due to the complexity of the task, make the assessment of precision and recall necessary. The binning process can be **taxonomy dependent**, relying on a database to aggregate the sequences, or **independent**. The first approach has many of the issues shared with the direct analysis of metagenomic read data as it’s highly dependent on the database use. 

The independent approach has the benefit of not relying on a database, instead using the composition of each sequence and coverage profiles to cluster together sequences that might belong to the same organism. These algorithms don’t require prior knowledge about the genomes in a given sample, relying instead on features inherent to the sequences in the sample. Although most binning softwares can work with single metagenomic samples, most make use of differential coverage of multiple samples to improve the binning process (Sedlar et al., 2016). This relies on the underlying principle that different organisms will be present in different proportions, therefore having a different coverage profile in the sample, which is not always true. Enrichment steps prior to sequencing, such as probes, can significantly alter the abundance profile making it more homogeneous, which might have profound consequences in the binning performance. The sampling as sequencing also introduce biases in the relative abundance of the organisms represented, either masking low abundance ones or over-representing certain taxas. 

The binning process represents a huge constraint in the analysis of metagenomic data. It’s an essential step to be able to retrieve information with a taxonomic context but thresholds need to be met on the purity and completeness of the bins obtained. Taxonomic independent are usually preferred as the results obtained aren’t constrained by the database used, but regardless of the methodology used, the results obtained are often unreliable and unreproducible due to the low recall of these methods (Sczyrba et al., 2017). The latest binning comparison happened in 2017, with the first CAMI challenge using simulated data sets, generated from ∼700 newly sequenced microorganisms and ∼600 novel viruses and plasmids and representing common experimental setups, where maxbin2 was considered superior (Sczyrba et al, 2017). To the best of our knowledge, no assessment of the influence of complexity of the communities nor the variability of proportion within the communities was performed. 


## Methods

#### Simulated dataset
With the [M3S3 tool](http://medweb.bgu.ac.il/m3s3/), a simulation sample was obtained made up of the Zymobiomics 
microbial community standasrd species of bacteria and yeast. The sample has the following composition (obtained with 
Kraken2 with the minimkaken2_V1 database):

#### Zymos community Standard
Two commercially available mock communities containing 10 microbial species (ZymoBIOMICS Microbial Community Standards) 
were sequences by [Nicholls et al. 2019](https://academic.oup.com/gigascience/article/8/5/giz043/5486468). Shotgun 
sequencing of the Even and Log communities was performed with the same protocol, with the exception that the Log 
community was sequenced individually on 2 flowcell lanes and the Even community was instead sequenced on an Illumina 
MiSeq using 2×151 bp (paired-end) sequencing. They are available under accession numbers 
[ERR2984773](https://www.ebi.ac.uk/ena/data/view/ERR2984773) (even) and 
[ERR2935805](https://www.ebi.ac.uk/ena/data/view/ERR2935805).

#### Real dataset
As a real metagenomic dataset, the three CAMI synthetic metagenome challenge sets was used (Low (L), Medium with 
differential abundance (M), and High complexity with time series (H)). These datasets are publicly available in 
[GigaDB](http://gigadb.org/dataset/100344). 


### Binning tools and Commands

### Assessing Binning success

## Results

## References

* Angers-Loustau, A., Petrillo, M., Bengtsson-Palme, J., Berendonk, T., Blais, B., Chan, K.-G., et al. (2018). The challenges of designing a benchmark strategy for bioinformatics pipelines in the identification of antimicrobial resistance determinants using next generation sequencing technologies. F1000Research 7, 459. doi:10.12688/f1000research.14509.1.
Benson, D. A. (2004). GenBank. Nucleic Acids Research 33. doi:10.1093/nar/gki063.
* Couto, N. et al. Critical steps in clinical shotgun metagenomics for the concomitant detection and typing of microbial pathogens. Sci. Rep. 8, 13767 (2018).
* Gruening, B., Sallou, O., Moreno, P., Leprevost, F. D. V., Ménager, H., Søndergaard, D., et al. (2018). Recommendations for the packaging and containerizing of bioinformatics software. F1000Research 7, 742. doi:10.12688/f1000research.15140.1.
Jia, B., Raphenya, A. R., Alcock, B., Waglechner, N., Guo, P., Tsang, K. K., et al. (2016). CARD 2017: expansion and model-centric curation of the comprehensive antibiotic resistance database. Nucleic Acids Research 45. doi:10.1093/nar/gkw1004.
* Sczyrba, A., Hofmann, P., Belmann, P., Koslicki, D., Janssen, S., Dröge, J., et al. (2017). Critical Assessment of Metagenome Interpretation—a benchmark of metagenomics software. Nature Methods 14, 1063–1071. doi:10.1038/nmeth.4458.
* Sedlar, K., Kupkova, K., and Provaznik, I. (2017). Bioinformatics strategies for taxonomy independent binning and visualization of sequences in shotgun metagenomics. Computational and Structural Biotechnology Journal 15, 48–55. doi:10.1016/j.csbj.2016.11.005.


## Authors

* Inês Mendes 
    * Instituto de Microbiologia, Instituto de Medicina Molecular, Faculdade de Medicina, Universidade de Lisboa, 
    Lisboa, Portugal; 
    * University of Groningen, University Medical Center Groningen, Department of Medical Microbiology and Infection 
    Prevention, Groningen, The Netherlands
