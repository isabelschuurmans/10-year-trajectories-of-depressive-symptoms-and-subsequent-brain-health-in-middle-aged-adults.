module load R/3.6.3
module load freesurfer/6.0.0-Mijnwerkplek

R
library(QDECR)
setwd("/home/r042647/depression_mri")
vw <- readRDS('rh.depression_clus.M2.thickness/rh.depression_clus.M2.thickness.rds')
freeview(vw, stack = 'clusB')

# load:
# - rh.pial as surface
# binary curvature
# overlay stack 2 ... masked
#---- if R doesnt work
freeview --surface /mnt/data/genr/mrdata/GenR_parent/bids/derivatives/freesurfer/6.0.0/fsaverage/surf/rh.inflated:overlay=/dev/shm/rh.depression_clus.M2.thickness.stack2.temp_mgh.mgh:overlay_method=linearopaque:overlay_threshold=0.0545253418385983,0.0791986957192421
