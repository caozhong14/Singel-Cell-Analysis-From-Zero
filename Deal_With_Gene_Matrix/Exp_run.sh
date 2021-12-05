#!/bin/bash
# ./Exp_Prepare_InstallPackages.R # if failed, please install packages manually.

# ./Exp_LoadAndBatch.r # Skip or run for differen intergration comparision

# ./Exp_Cluster.r  | tee logs/log_Exp_Cluster.txt # Set False in line 173

# ./plot_Figure1.R -i "../output/ASSHD_pbmc_orig.rds"  \
#                    -c "integrated_snn_res.0.8" \
#                    -o "../Figure1_Cluster_Primary8/" | tee logs/log_plotFigure1_Cluster_Primary8.txt 


# ./plot_Figure1.R -i "../output/ASSHD_pbmc_orig.rds"  \
#                    -c "integrated_snn_res.0.4" \
#                    -o "../Figure1_Cluster_Primary4/" | tee logs/log_plotFigure1_Cluster_Primary4.txt 

# ./plot_Figure1.R -i "../output/ASSHD_pbmc_orig.rds"  \
#                    -c "integrated_snn_res.0.6" \
#                    -o "../Figure1_Cluster_Primary6/" | tee logs/log_plotFigure1_Cluster_Primary6.txt 

######################################################################################
## If lucky, you will identify cell type from these Figure1_Cluster_Primary4/6/8    ##
## If not, you can analyze the features further more, codes are prepared            ##
## Codes are in Exp_RenameIdents.rmd, please run it and analyze results             ##
## If some clusters has too much markers, Refer to Quality Control in Exp_Cluster.r ##
######################################################################################

# ./Exp_Cluster.r # Set True in line 173 to identify cell types

# ./plot_Figure1.R -i "../output/ASSHD_pbmc_cluster.rds"  \
#                    -c "celltypeI" \
#                    -o "../Figure1_Cluster_I/" | tee logs/log_plotFigure1_Cluster_I.txt 

# ./plot_Figure1.R -i "../output/ASSHD_pbmc_subcluster_B.rds"  \
#                    -c "celltypeII_B" \
#                    -o "../Figure1_Cluster_II_B/" | tee logs/log_plotFigure1_Cluster_II_B.txt 

# ./plot_Figure1.R -i "../output/ASSHD_pbmc_subcluster_Bcom.rds"  \
#                    -c "celltypeII_B" \
#                    -o "../Figure1_Cluster_II_Bcom/" | tee logs/log_plotFigure1_Cluster_II_Bcom.txt 

# ./plot_Figure1.R -i "../output/ASSHD_pbmc_subcluster_T.rds"  \
#                    -c "celltypeII_T" \
#                    -o "../Figure1_Cluster_II_T/" | tee logs/log_plotFigure1_Cluster_II_T.txt 

./plot_Figure1.R -i "../output/ASSHD_pbmc_subcluster_Tcom.rds"  \
                   -c "celltypeII_T" \
                   -o "../Figure1_Cluster_II_Tcom/" | tee logs/log_plotFigure1_Cluster_II_Tcom.txt 

# ./plot_Figure1.R -i "../output/ASSHD_pbmc_orig.rds"  \
#                    -c "integrated_snn_res.0.8" \
#                    -o "../Figure1_Cluster_II_I_res0.8/" | tee logs/log_plotFigure1_Cluster_I_res0.8.txt 

# ./plot_Figure1.R -i "../output/ASSHD_pbmc_orig.rds"  \
#                    -c "integrated_snn_res.1" \
#                    -o "../Figure1_Cluster_II_I_res1.0/" | tee logs/log_plotFigure1_Cluster_I_res1.0.txt 


# ./script_go_kegg.R -i "../output/ASSHD_pbmc_cluster_com.rds" -l 0.3 -o "../Figure2_DEG_0.3" | tee logs/log_ASSHD_cluster_I_go_kegg_fc.3.txt 

# ./script_go_kegg.R -i "../output/ASSHD_pbmc_cluster_com.rds" -l 0.5 -o "../Figure2_DEG_0.5" | tee logs/log_ASSHD_cluster_I_go_kegg_fc.5.txt 

# ./script_go_kegg.R -i "../output/ASSHD_pbmc_cluster_com.rds" -l 0.8 -o "../Figure2_DEG_0.8" | tee logs/log_ASSHD_cluster_I_go_kegg_fc.8.txt 
