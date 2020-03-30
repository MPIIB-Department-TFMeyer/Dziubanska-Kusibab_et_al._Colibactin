#setwd("./Code/3_DSB-shape-characteristics/Multivariate/")
load("../../../Data/Shape_correlations/Multivariate/6nt_data_with_stiffness_and_PCA.Rdata")
library(rgl)

inp = as.matrix(d_param_full_6nt[, 31:108])
inp = inp[order(-d_param_full_6nt$cl),]
rownames(inp) = d_param_full_6nt$Nt

other_col = rgb(0.9, 0.9, 0.9, 0.1)


sel_HN_pattern = "AATTTT|AAATTT|AAAATT|AATATT"

# Selected HN + enrichment in central HN highlighted
col = with(d_param_full_6nt, ifelse(grepl(sel_HN_pattern,Nt),ifelse(cl>0.2,"red", "green"),ifelse(cl>0.2,"blue",other_col)))

col = col[order(-d_param_full_6nt$cl)]

plot3d(pca$x[,1], pca$x[,2], pca$x[, 3], col=col, type="n", xlab="PC 1", ylab="PC 2", zlab="PC 3")

# Define a maximum number of data points per volume in the 3D-plot. For voxels with more points sample
# randomly from the existing points in the voxels. 
dim_ranges = apply(pca$x[,1:3], 2, range)
pos_discrete = pca$x[,1:3]
for (cc in colnames(dim_ranges)) {
  pos_discrete[, cc] = cut(pos_discrete[, cc], breaks=seq(dim_ranges[1,cc], dim_ranges[2,cc], by=20/diff(dim_ranges[,1])))
}
pos_volume_id = apply(pos_discrete, 1, paste, collapse="_")
per_volume_id_counts = table(pos_volume_id)
max_points_per_voxel = 50
high_count_per_volume_id_voxels = sum(per_volume_id_counts>max_points_per_voxel)

selected_points = tapply(names(pos_volume_id), factor(pos_volume_id), 
                         function(x) {if(length(x)>max_points_per_voxel) sample(x,max_points_per_voxel) else x},
                         simplify=F)

#subsample = sample.int(nrow(d_param_full_6nt), 50000)
#subsample = 1:nrow(d_param_full_6nt)
subsample = (d_param_full_6nt$Nt %in% unlist(selected_points))
tmp = d_param_full_6nt[subsample,]
size=with(tmp, ifelse(grepl(sel_HN_pattern,Nt),2,0.2))/2
size = size[order(-d_param_full_6nt$cl)]
inc_flag = ifelse(size==max(size), T, F)

rgl.spheres(pca$x[subsample,1], pca$x[subsample,2], pca$x[subsample, 3], r = 0.35, color = col[subsample])
view3d(-75, -65, zoom=0.5, fov=135)

rgl.snapshot(filename = "../../../Results/3_DSB-shape-characteristics/Snapshot_3D_v2020_01-28.png",top = TRUE)
