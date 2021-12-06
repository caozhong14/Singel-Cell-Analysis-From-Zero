## Workspace and Data

```
.
├── data -> ../scRNAdataset
├── Exps
├── scRNA-tools
└── Singel-Cell-Analysis-From-Zero
    ├── Easy_Way_To_Use_Cell_Ranger
    ├── LICENSE
    ├── README.html
    └── README.md
```

## Step 1. Install scRNA-tools
Distributor ID:	Ubuntu
Description:	Ubuntu 16.04.7 LTS
Release:	16.04
Codename:	xenial

```
cd scRNA-tools
```

|  Software   | version  | function | link
|  ----  | ----  | ----  | ----  |
| FastQC  | v0.11.9 | Quality Control | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| Cell Ranger  | 6.0.0 | Cell Counter | https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger |


### Install FastQC_v0.11.9.zip

Install Jave
```Terminal
sudo apt-get update
sudo apt-get install default-jre
```

Download and install fastqc.
```Terminal
cd /YOUR/WORKSPACE/PATH/
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
cd FastQC/
chmod 755 fastqc
sudo ln -s /YOUR/WORKSPACE/PATH/FastQC/fastqc /usr/local/bin/fastqc
```

### Install Cell Ranger - 6.0.0
Download Cell Ranger
```curl 
curl -o cellranger-6.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.0.0.tar.gz?Expires=1617326791&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MTczMjY3OTF9fX1dfQ__&Signature=F47OynLO8C-uo4OOYegoHxTu4bsML0OhgKJnvV1lPbNFt6OWSix13wnhfhWEZ0MbN9w8DPCSJKqO5fLBxkch8bkLNFpYi6lVpSD278QP7MtUwDIb6rU6ijJXdwH-c2LILA9E1te3mWSeUgoZtk10D4LRvdSSgifTYdd96GLHDANkdSIjv3uT5l53-zhidSSlR7D9cuGU9Da706NtKek5v0m-FWDiI25Ngof0MMjO~PQbYx-cAxVS8Me8BexG9TvOQMLYCPs5qSGl5vujVSHEye~freZ1UJvJOBV~9LacVLPgX~2Q2aDHdAg6ZcOcKvn9H~tm0YBhGa5CdHdXU8YJsA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
```

```wget 
wget -O cellranger-6.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.0.0.tar.gz?Expires=1617376081&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MTczNzYwODF9fX1dfQ__&Signature=kcU-bYLzAEqN7gjyuomIpkXIRT1fBvWHYWiodTXCu3SJ44ajhgpU0VMXkOrZk1X6oJZPu6OuwsxcEyzE5fVJCWcwwTLrTAY6gW~-rqM5YLJqFiNKaiurSoHsDI-w8NqIyoBOrc1EdGIwdSUgeM5WBGEFpzWGaOCsaqAGCBEHeGH3GGCm7uqCG8ZpnZWC6oKKJ~wyVl9D~j8cs6xiDSRsA4u2Juy3lSA5h44PzReWoi2BaehTBLYqxp6AHyCcqsvB5RxldDrdwiMrgmygERO2bz6xUIz81zTjM9HWJsE6PO7HPrIZJJtZAZFU2FaS0y0AKjdhp-Q8E1cDbUMAoqW3AQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
```

```
tar -zxvf cellranger-6.0.0.tar.gz 
```

Download Human reference (GRCh38)
```curl 
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
```

```wget 
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
```

```
tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz
```

### Install Loupe Cell Browser 5.0.1

Official Website.
https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest?

### Install R & Seurat

```
wget https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/base/R-4/R-4.1.2.tar.gz
tar -zxvf R-4.1.2.tar.gz
```
``` R install
cd R-4.1.2
./configure --prefix=/path/R/4.1.2 
sudo apt-get install libcurl4-openssl-dev
make && make install
```
Or install with apt-get command.
```
sudo apt-get install --no-install-recommends r-base
```

## Got Filtered Feature Matrix with Cell Ranger workflow

The main processes of Cell Ranger are mainly: mkfastq, count, aggr, and reanalyze. There are also some small tools, such as mkref, mkgtf, upload, sitecheck, mat2csv, vdj, mkvdjref, testrun. Please refer to official documents.

https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ov

### Preparetion for environment

Set your workspace as /home/caozhong/trunk/SingleCellAnalysis. Including data, codes, tools and documents.

```
.
├── data -> ../scRNAdataset
├── Exps
├── scRNA-tools
└── Singel-Cell-Analysis-From-Zero
    ├── Easy_Way_To_Use_Cell_Ranger
    ├── LICENSE
    ├── README.html
    └── README.md
```

``` 
export SCRNA_HOME='/home/caozhong/trunk/SingleCellAnalysis'
alias cellranger=cellranger
```

```
cellranger sitecheck > sitecheck.txt
vim sitecheck.txt
cellranger testrun --id=check_install
```
In order to separate data, tools, tutorials and experimental results, we recommend that each experiment be conducted in Exps to avoid changes to data, tools, and tutorials.

```
cd Exps
```


### mkfastq 

We don't use mkfastq, because the fastq data is downloaded from the website or the company. They usually gives fastq data. We introduce mkfastq here in case you need to prepare your own data:
``` bash
cellranger mkfastq --id=bcl \
                   --run=/path/to/bcl \
                   --samplesheet=samplesheet-1.2.0.csv \
                   --fastqs=/PATH/TO/PROJECT_FOLDER \ # ./201292A_T1PBMC_5
                   --sample=201292A_T1PBMC_5

# id: output name
# run: input BCL File Folder

```
When the BCL filename is not standardized, we need to rename like, rename 201292D_T3PBMC_5_2_1.R2.fq.gz as 201292D_T3PBMC_5_S1_L002.R2_001.fastq.gz. We can
```
rename "s/C_5_/C_5_S1_L00/" *
rename "s/1_R1.fq.gz/R1_001.fq.gz/" *
rename "s/1_R2.fq.gz/R2_001.fq.gz/" *
rename "s/fq/fastq/"
```


### cellranger 
For cell count. There are many internal processes: comparison, quality control, quantification and so on. Just use it simply.

``` bash 
cellranger count --id=tmp_T1PBMC_5_S1 \

--transcriptome=/home/caozhong/mytrunk/SingleCellAnalysis/scRNA-tools/refdata-gex-GRCh38-2020-A \

--fastqs="/home/caozhong/mytrunk/SingleCellAnalysis/data/T1PBMC_202104_data/KYSY-2764-JD-YX-2020-1292-JSFU-04/origData/201292A_T1PBMC_5" \

--sample=201292A_T1PBMC_5 \

--expect-cells=5000 \

--nosecondary
```

Finally, we got data_filtered_feature_bc_matrix and move data to Singel-Cell-Analysis-From-Zero/data_filtered_feature_bc_matrix.

```

├── all.md5.saved
├── ASS5HD3
│   ├── ASS-ILD1
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── ASS-ILD2
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── ASS-ILD3
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── ASS-ILD4
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── ASS-ILD5
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── HD1
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── HD2
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   └── HD3
│       ├── barcodes.tsv.gz
│       ├── features.tsv.gz
│       └── matrix.mtx.gz

```