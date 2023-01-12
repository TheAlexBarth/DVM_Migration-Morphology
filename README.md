# DVM Morphology and Seasonality

## R

All rscripts will produce summary data which are stored in the /data directory.
All files are stored as RDS which provide the smallest single file structure for hierarchical data. However, there are some challenges version control on github. Files about 100mb are prohibited. So for some files I label as _large.rds
These files are placed on git large file storage. There is no big difference either way. To access raw data, contact Josh Stone (stone@sc.edu). To access formatted data in a non-R friendly format, contact Alex Barth (abarth1014@gmail.com).

Files and general workflow:
Several of these files rely on a data structure created in a prior file. All data are created but to run from the raw data, follow this sequence.

-   mgmt_00_modis-data-import: imports raw modis data and formats as a list
-   mgmt_00_proile-org: imports raw data and filters/formats as needed
-   mgmt_01_geog-filter: filters ctd and uvp to matching gps locations; creates a gps location file
-   mgmt_ctd_interp-profiles: creates interpolation grid for plotting ctd profiles across the environmental gradients
-   mgmt_02_PCA-copepods: runs pca and k-means on copepods. Produces pca results (data/02_cope-pca-res) and cluster groups (data/02_cluster-copepods)
-   mgmt_03_cluster_conc: calculates cluster concentration for each cast (data/03_cluster_conc)
-   mgmt_03_wmd-boot-all: generates resampled distributions for wmd of both all data (data/03_wmd-all-boot)
-   mgmt_04_occ-model: runs occurance model with all plankton data to visualize volume sampled effects creates (data/04_full-col-model-*.rds)
-   mgmt_05_cluster-occ-model: runs occurence model on individual clusters creates (data/05_*.rds)
-   mgmt_06_environ-occ-model: runs environmental occurence model and creates (data/06_*.rds)

-   jags_04_occ-detection-only: a jags script to run the detection only model (runs on 5 too)
-   jags_06_occ-model: jags script to run environmental data occurence model