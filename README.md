# variant2literature (v2l)

Extract and normalize variants from academic papers in xml, pdf, doc, docx, xlsx, csv formats.

## Required Data and Packages
##### CRF++:
- download <a href=https://taku910.github.io/crfpp/>CRF++.0.58.tar.gz</a> \
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
- https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
- https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz

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
### Run variant2literature in Docker
#### Prerequisites
- Linux OS
- Docker 18.09.0 or higher
- CUDA 8.0 or higher
- nvidia-docker  

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


### Directly run in Linux
#### Setup

```sh
apt-get install -y software-properties-common
add-apt-repository ppa:deadsnakes/ppa
apt-get update
apt-get install -y \
        build-essential cmake \
        python3.6-dev python3-pip python3-tk \
        libpoppler-cpp-dev libmagic-dev libxrender-dev \
        libsm6 libxext6 libglib2.0-0 \
        libreoffice poppler-utils
```
##### Install python and required package
```sh
ln -s /usr/bin/python3.6 /usr/local/bin/python
python -m pip install -U pip==18.1
pip install torch==0.4.1
# If you have CUDA 9.2, please use the following command to install pytorch instead
# pip install http://download.pytorch.org/whl/cu92/torch-0.4.1-cp36-cp36m-linux_x86_64.whl

pip install -r requirements.txt \
    && python -c "import nltk; nltk.download('punkt'); nltk.download('averaged_perceptron_tagger')" 
```
##### Initialize CRF++
```sh
cp CRF++-0.58.tar.gz /opt/
cd /opt && tar zxvf CRF++-0.58.tar.gz
cd /opt/CRF++-0.58 \
    && ./configure && make && make install && cd python \
    && cp /opt/CRF++-0.58/crfpp.h . \
    && python setup.py build && ldconfig \
    && python setup.py install
```
##### Initialize table detector
```sh
cd table_detector/lib && bash make.sh
```
##### Install mysql
```sh
apt-get install mariadb-server
service mysql start
```

##### Change mysql root password
```sh
mysql_secure_installation
# then enter default password `s8fjYJd92oP`
```
If you get en error like `1698, "Access denied for user 'root'@'localhost'"`, please set the root user to use the mysql_native_password plugin.
```
mysql> USE mysql;
mysql> UPDATE user SET plugin='mysql_native_password' WHERE User='root';
mysql> FLUSH PRIVILEGES;
mysql> exit;
```
then restart mysql
```sh
service mysql restart
```
##### Load mysql data
```sh
ln -s ./ /app

export MYSQL_HOST=127.0.0.1
export MYSQL_PORT=3306
export MYSQL_ROOT_PASSWORD=s8fjYJd92oP

cd mysqldb
python models.py
```
##### Run table detector
```sh
export CUDA_VISIBLE_DEVICES=0
export NUM_TABLE_DETECTORS=1
export LOAD_BALANCER_HOST='localhost'

cd table_detector && python table_detector.py
```

##### Index papers
Put paper directories in `input/`, then execute
```sh
python main.py --n-process 1 --input input/
```
If your input files are plain text, or you're running on a device without GPU, please add `--no-table-detect` to disable the table detector.   
The results will be saved in mysql database, please use `query.py` to query or use SQL command directly. For example:
```
mysql> USE gene;
mysql> SELECT * FROM var_pmid WHERE _id='<paper_directory_name>';
```

##### Query
```sh
python query.py
```

## License

This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

The [fasterRCNN implementation](https://github.com/jwyang/faster-rcnn.pytorch) used here is written by [Jianwei Yang](https://github.com/jwyang) and [Jiasen Lu](https://github.com/jiasenlu).
