# variant2literature

Extract and normalize variants from academic papers in xml, pdf, doc, docx, xlsx, csv formats.

## Prerequisites
- Linux OS
- Docker 18.09.0 or higher
- CUDA 8.0 or higher
- nvidia-docker

## Required Data and Packages

##### CRF++:
- download <a href=https://taku910.github.io/crfpp/>CRF++.0.58.tar.gz<a> \
put CRF++.0.58.tar.gz in `variant2literature/`


#### download following files and put them in `variant2literature/models/`

##### FasterRCNN model:
- https://www.dropbox.com/s/g980k8hpqj1q8cn/faster_rcnn.pth

##### UCSC tables (hg19):
- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/ncbiRefSeq.txt.gz
- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/ensGene.txt.gz
- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/knownGene.txt.gz
- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/kgAlias.txt.gz
- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/ensemblToGeneName.txt.gz
- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp150.txt.gz

##### NCBI gene_info
- ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
- ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz

##### ucsc.hg19.fasta
- download from ucsc and convert it to fasta format or download from GATK bundle and decompress it. \
rename it to `ucsc.hg19.fasta` (if the filename is not `ucsc.hg19.fasta`).

##### tmVar
- download <a href="https://www.ncbi.nlm.nih.gov/research/bionlp/Tools/tmvar/">tmVar 2.0</a> \
copy `tmVarJava/CRF/MentionExtractionUB.Model` to `variant2literature/models/`

##### GNormPlus
- download <a href=https://www.ncbi.nlm.nih.gov/research/bionlp/Tools/gnormplus/>GNormPlus</a> \
copy `GNormPlusJava/Dictionary/GNR.Model` to `variant2literature/models/` and \
copy `GNormPlusJava/Dictionary/PT_CTDGene.txt` to `variant2literature/models/`

## Usage

#### Setup

- build docker image by `make build`
- compile fasterRCNN by `make compile`
- start docker container by `make run`
- start mysql docker container by `make run-db`
- load data into database by `make load-db` (run only once unless MYSQL_VOLUME is changed)

#### Index Papers

- put paper directories in input/
- run `make index`
- query by `make query` or `make query OUTPUT_FILE=output.txt`

#### Delete Indexes

- run `make truncate`

#### Stop and Remove Docker Container

- run `make rm`
- run `make rm-db`

## License

This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details.
