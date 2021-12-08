# Single-Cell-Analysis-From-Zero

Codes for paper "Single-cell transcriptomics reveals peripheral immune responses in anti-synthetase syndrome-associated interstitial lung disease".

## Got Filtered Feature Matrix with Cell Ranger workflow

Please refer to the file "Easy_Way_To_Use_Cell_Ranger/Cell Ranger workflow.md" to get matrices.

Patient: ASS-ILD1, ASS-ILD2, ASS-ILD3, ASS-ILD4, ASS-ILD5
Control: HD1, HD2, HD3

We have uploaded our raw data and processed matrices to GEO database geo@ncbi.nlm.nih.gov, expect to receive an email from GEO curators with GEO accession numbers within 5 business days. And the alternative download link is https://www.aliyundrive.com/s/ECGrvFb9DeK. Please contact caozhong@mails.tsinghua.edu.cn to get the download code.

We also prepared the matrices, you can directly download with the following commands. 

```bash
wget http://http://82.157.178.197/f/07240b4a0cba4c77adfb/?dl=1 -O ASS5HD3.zip
unzip ASS5HD3.zip
mv ASS5HD3 data_filtered_feature_bc_matrix/
```

## Deal with Gene Matrix

Please refer to the commands "Deal_With_Gene_Matrix/Exp_run.sh".


###  Introduce to Script

1. Exp_Prepare_InstallPackages.r 

- Install R packages

2. Exp_LoadAndBatch.r

- Run it to observe whether we should do batch correction.

3. Exp_Cluster.r

- The basic single cell analysis workflow code.

4. Exp_RenameIdents_Cluster_I.rmd/II_B.rmd/II_T.rmd

- Tips for cluster and gene analysis.

5. plot_Figure1.R

- Draw figures for cluster.

6. script_go_kegg.R

- The basic GO and KEGG analysis code.

### Command to use scripts

You can run the commands in Deal_With_Gene_Matrix/Exp_run.sh one-by-one or directly run the Exp_run.sh script.

```bash
cd Deal_With_Gene_Matrix
bash Exp_run.sh
```


Then you will get results saved in the folder Singel-Cell-Analysis-From-Zero/output.
And you can get figures in the folders Singel-Cell-Analysis-From-Zero/Figure*

For more detailed analysis, please refer to our paper: "Single-cell transcriptomics reveals peripheral immune responses in anti-synthetase syndrome-associated interstitial lung disease".

If there are any questions about the data, the codes or the results, please do not hesitate to contact us caozhong14@mails.tsinghua.edu.cn.












