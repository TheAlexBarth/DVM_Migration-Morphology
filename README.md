# DVM Morphology and Seasonality

## R

All rscripts will produce summary data which are stored in the /data directory.
All files are stored as RDS which provide the smallest single file structure for hierarchical data. However, there are some challenges version control on github. Files about 100mb are prohibited. So for some files I label as _large.rds
These files are placed on git large file storage. There is no big difference either way. To access raw data, contact Josh Stone (stone@sc.edu). To access formatted data in a non-R friendly format, contact Alex Barth (abarth1014@gmail.com).

Files and general workflow:

-   mgmt_00_modis-data-import: imports raw modis data and formats as a list
-   mgmt_00_proile-org: imports raw data and filters/formats as needed
-   mgmt_01_geog-filter: filters ctd and uvp to matching gps locations; creates a gps location file



-   mgmt_ctd_02_interp-profiles: creates interpolated ctd data (data/ctd_02_interp-grid). note this would typically be a plot\_ but can't store large rasters. So plots are coded in quarto_dev
-   mgmt_02_PCA-copepods: runs pca and k-means on copepods. Produces pca results (data/02_cope-pca-res) and cluster groups (data/02_cluster-copepods)
-   mgmt_03_cluster_conc: calculates cluster concentration for each cast (data/03_cluster_conc)
-   mgmt_04_wmd-boot-all: generates resampled distributions for wmd of both all data (data/04_wmd-all-boot) and per cruise (data/04_wmd-by-cruise.rds)