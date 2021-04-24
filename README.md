# Single Cell RNA-Seq Analysis of Pancreatic Cells

The pancreas is a complex organ comprised of a diverse set of cell types. Proper function of the pancreas is required to maintain healthy metabolism, and pancreatic dysfunction leads to serious illnesses. Baron et al performed single cell RNA sequencing in a set of post-mortem human donor pancreatic cells from four subjects and two mouse models to better understand the cellular diversity in the pancreas. Analysis of the data identified previously known cell types as well as rare and novel cell type subpopulations, and created a more detailed characterization of the diversity of those cell types. In this project, we will attempt to replicate their primary findings using current analytical methodology and software packages.

Project Goals:
* Process the barcode reads of a single cell sequencing dataset
* Perform cell-by-gene quantification of UMI counts
* Perform quality control on a UMI counts matrix
* Analyze the UMI counts to identify clusters and marker genes for distinct cell type populations
* Ascribe biological meaning to the clustered cell types and identify novel marker genes associated with them

Original Analysis:
Baron, Maayan, Adrian Veres, Samuel L. Wolock, Aubrey L. Faust, Renaud Gaujoux, Amedeo Vetere, Jennifer Hyoje Ryu, et al. 2016. “A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-Cell Population Structure.” Cell Systems 3 (4): 346–60.e4. PMID: 27667365

## Contributors

* Daisy Wenyan Han daisyhan@bu.edu
* Divya Sundaresan divyas3@bu.edu
* Alec Jacobsen aggjacob@bu.edu
* Emmanuel Saake esaake@bu.edu

## Repository Contents and Suggested Workflow

* `data_curator_qsubs` - Using the raw reads, in FASTQ file format, a whitelist of barcodes is generated to properly filter poor-quality samples by UMI count. A raw UMI count matrix is then generated using the Salmon Alevin program, available on the command line.
* `programmer.R` - The UMI count matrix generated above is loaded, and processed using the Seurat standard pre-processing workflow, as detailed by the Satija lab. Low quality reads are filtered, and cells are clustered into subpopulations.
* `analyst.R` - The cell subpopulations are classified into distinct cell types, using the marker genes provided in the Supplementary Data of Baron et al.. Marker genes for each of these cell types are retained, with a list of novel marker genes exported for further analysis.
* `Biologist.R` - The novel marker genes identified in the step above are examined for gene set enrichment, identifying key pathways that are over- or under-expressed, by cell type.
